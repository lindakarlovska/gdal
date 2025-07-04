/******************************************************************************
 *
 * Project:  Mapinfo Image Warper
 * Purpose:  Implementation of one or more GDALTrasformerFunc types, including
 *           the GenImgProj (general image reprojector) transformer.
 * Author:   Frank Warmerdam, warmerdam@pobox.com
 *
 ******************************************************************************
 * Copyright (c) 2002, i3 - information integration and imaging
 *                          Fort Collin, CO
 * Copyright (c) 2008-2013, Even Rouault <even dot rouault at spatialys.com>
 * Copyright (c) 2021, CLS
 *
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#include "cpl_port.h"
#include "gdal_alg.h"
#include "gdal_alg_priv.h"

#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <limits>
#include <utility>

#include "cpl_conv.h"
#include "cpl_error.h"
#include "cpl_list.h"
#include "cpl_minixml.h"
#include "cpl_multiproc.h"
#include "cpl_string.h"
#include "cpl_vsi.h"
#include "gdal.h"
#include "gdal_priv.h"
#include "ogr_core.h"
#include "ogr_spatialref.h"
#include "ogr_srs_api.h"

CPL_C_START
void *GDALDeserializeGCPTransformer(CPLXMLNode *psTree);
void *GDALDeserializeTPSTransformer(CPLXMLNode *psTree);
void *GDALDeserializeGeoLocTransformer(CPLXMLNode *psTree);
void *GDALDeserializeRPCTransformer(CPLXMLNode *psTree);
void *GDALDeserializeHomographyTransformer(CPLXMLNode *psTree);
CPL_C_END

static CPLXMLNode *GDALSerializeReprojectionTransformer(void *pTransformArg);
static void *GDALDeserializeReprojectionTransformer(CPLXMLNode *psTree);

static CPLXMLNode *GDALSerializeGenImgProjTransformer(void *pTransformArg);
static void *GDALDeserializeGenImgProjTransformer(CPLXMLNode *psTree);

static void *GDALCreateApproxTransformer2(GDALTransformerFunc pfnRawTransformer,
                                          void *pRawTransformerArg,
                                          double dfMaxErrorForward,
                                          double dfMaxErrorReverse);

/************************************************************************/
/*                            GDALIsTransformer()                       */
/************************************************************************/

bool GDALIsTransformer(void *hTransformerArg, const char *pszClassName)
{
    if (!hTransformerArg)
        return false;
    // All transformers should have a GDALTransformerInfo member as their first members
    GDALTransformerInfo *psInfo =
        static_cast<GDALTransformerInfo *>(hTransformerArg);
    return memcmp(psInfo->abySignature, GDAL_GTI2_SIGNATURE,
                  strlen(GDAL_GTI2_SIGNATURE)) == 0 &&
           strcmp(psInfo->pszClassName, pszClassName) == 0;
}

/************************************************************************/
/*                          GDALTransformFunc                           */
/*                                                                      */
/*      Documentation for GDALTransformFunc typedef.                    */
/************************************************************************/

/*!

\typedef typedef int (*GDALTransformerFunc)( void *pTransformerArg, int
bDstToSrc, int nPointCount, double *x, double *y, double *z, int *panSuccess );

Generic signature for spatial point transformers.

This function signature is used for a variety of functions that accept
passed in functions used to transform point locations between two coordinate
spaces.

The GDALCreateGenImgProjTransformer(), GDALCreateReprojectionTransformerEx(),
GDALCreateGCPTransformer() and GDALCreateApproxTransformer() functions can
be used to prepare argument data for some built-in transformers.  As well,
applications can implement their own transformers to the following signature.

\code
typedef int
(*GDALTransformerFunc)( void *pTransformerArg,
                        int bDstToSrc, int nPointCount,
                        double *x, double *y, double *z, int *panSuccess );
\endcode

@param pTransformerArg application supplied callback data used by the
transformer.

@param bDstToSrc if TRUE the transformation will be from the destination
coordinate space to the source coordinate system, otherwise the transformation
will be from the source coordinate system to the destination coordinate system.

@param nPointCount number of points in the x, y and z arrays.

@param[in,out] x input X coordinates.  Results returned in same array.

@param[in,out] y input Y coordinates.  Results returned in same array.

@param[in,out] z input Z coordinates.  Results returned in same array.

@param[out] panSuccess array of ints in which success (TRUE) or failure (FALSE)
flags are returned for the translation of each point. Must not be NULL.

@return TRUE if all points have been successfully transformed (changed in 3.11,
previously was TRUE if some points have been successfully transformed)

*/

/************************************************************************/
/*                      GDALSuggestedWarpOutput()                       */
/************************************************************************/

/**
 * Suggest output file size.
 *
 * This function is used to suggest the size, and georeferenced extents
 * appropriate given the indicated transformation and input file.  It walks
 * the edges of the input file (approximately 20 sample points along each
 * edge) transforming into output coordinates in order to get an extents box.
 *
 * Then a resolution is computed with the intent that the length of the
 * distance from the top left corner of the output imagery to the bottom right
 * corner would represent the same number of pixels as in the source image.
 * Note that if the image is somewhat rotated the diagonal taken isn't of the
 * whole output bounding rectangle, but instead of the locations where the
 * top/left and bottom/right corners transform.  The output pixel size is
 * always square.  This is intended to approximately preserve the resolution
 * of the input data in the output file.
 *
 * The values returned in padfGeoTransformOut, pnPixels and pnLines are
 * the suggested number of pixels and lines for the output file, and the
 * geotransform relating those pixels to the output georeferenced coordinates.
 *
 * The trickiest part of using the function is ensuring that the
 * transformer created is from source file pixel/line coordinates to
 * output file georeferenced coordinates.  This can be accomplished with
 * GDALCreateGenImgProjTransformer() by passing a NULL for the hDstDS.
 *
 * @param hSrcDS the input image (it is assumed the whole input image is
 * being transformed).
 * @param pfnTransformer the transformer function.
 * @param pTransformArg the callback data for the transformer function.
 * @param padfGeoTransformOut the array of six doubles in which the suggested
 * geotransform is returned.
 * @param pnPixels int in which the suggest pixel width of output is returned.
 * @param pnLines int in which the suggest pixel height of output is returned.
 *
 * @return CE_None if successful or CE_Failure otherwise.
 */

CPLErr CPL_STDCALL GDALSuggestedWarpOutput(GDALDatasetH hSrcDS,
                                           GDALTransformerFunc pfnTransformer,
                                           void *pTransformArg,
                                           double *padfGeoTransformOut,
                                           int *pnPixels, int *pnLines)

{
    VALIDATE_POINTER1(hSrcDS, "GDALSuggestedWarpOutput", CE_Failure);

    double adfExtent[4] = {};

    return GDALSuggestedWarpOutput2(hSrcDS, pfnTransformer, pTransformArg,
                                    padfGeoTransformOut, pnPixels, pnLines,
                                    adfExtent, 0);
}

static bool GDALSuggestedWarpOutput2_MustAdjustForRightBorder(
    GDALTransformerFunc pfnTransformer, void *pTransformArg, double *padfExtent,
    int /* nPixels*/, int nLines, double dfPixelSizeX, double dfPixelSizeY)
{
    double adfX[21] = {};
    double adfY[21] = {};

    const double dfMaxXOut = padfExtent[2];
    const double dfMaxYOut = padfExtent[3];

    // Take 20 steps.
    int nSamplePoints = 0;
    for (double dfRatio = 0.0; dfRatio <= 1.01; dfRatio += 0.05)
    {
        // Ensure we end exactly at the end.
        if (dfRatio > 0.99)
            dfRatio = 1.0;

        // Along right.
        adfX[nSamplePoints] = dfMaxXOut;
        adfY[nSamplePoints] = dfMaxYOut - dfPixelSizeY * dfRatio * nLines;
        nSamplePoints++;
    }
    double adfZ[21] = {};

    int abSuccess[21] = {};

    pfnTransformer(pTransformArg, TRUE, nSamplePoints, adfX, adfY, adfZ,
                   abSuccess);

    int abSuccess2[21] = {};

    pfnTransformer(pTransformArg, FALSE, nSamplePoints, adfX, adfY, adfZ,
                   abSuccess2);

    nSamplePoints = 0;
    int nBadCount = 0;
    for (double dfRatio = 0.0; dfRatio <= 1.01; dfRatio += 0.05)
    {
        const double expected_x = dfMaxXOut;
        const double expected_y = dfMaxYOut - dfPixelSizeY * dfRatio * nLines;
        if (!abSuccess[nSamplePoints] || !abSuccess2[nSamplePoints] ||
            fabs(adfX[nSamplePoints] - expected_x) > dfPixelSizeX ||
            fabs(adfY[nSamplePoints] - expected_y) > dfPixelSizeY)
        {
            nBadCount++;
        }
        nSamplePoints++;
    }

    return nBadCount == nSamplePoints;
}

static bool GDALSuggestedWarpOutput2_MustAdjustForBottomBorder(
    GDALTransformerFunc pfnTransformer, void *pTransformArg, double *padfExtent,
    int nPixels, int /* nLines */, double dfPixelSizeX, double dfPixelSizeY)
{
    double adfX[21] = {};
    double adfY[21] = {};

    const double dfMinXOut = padfExtent[0];
    const double dfMinYOut = padfExtent[1];

    // Take 20 steps.
    int nSamplePoints = 0;
    for (double dfRatio = 0.0; dfRatio <= 1.01; dfRatio += 0.05)
    {
        // Ensure we end exactly at the end.
        if (dfRatio > 0.99)
            dfRatio = 1.0;

        // Along right.
        adfX[nSamplePoints] = dfMinXOut + dfPixelSizeX * dfRatio * nPixels;
        adfY[nSamplePoints] = dfMinYOut;
        nSamplePoints++;
    }
    double adfZ[21] = {};

    int abSuccess[21] = {};

    pfnTransformer(pTransformArg, TRUE, nSamplePoints, adfX, adfY, adfZ,
                   abSuccess);

    int abSuccess2[21] = {};

    pfnTransformer(pTransformArg, FALSE, nSamplePoints, adfX, adfY, adfZ,
                   abSuccess2);

    nSamplePoints = 0;
    int nBadCount = 0;
    for (double dfRatio = 0.0; dfRatio <= 1.01; dfRatio += 0.05)
    {
        const double expected_x = dfMinXOut + dfPixelSizeX * dfRatio * nPixels;
        const double expected_y = dfMinYOut;
        if (!abSuccess[nSamplePoints] || !abSuccess2[nSamplePoints] ||
            fabs(adfX[nSamplePoints] - expected_x) > dfPixelSizeX ||
            fabs(adfY[nSamplePoints] - expected_y) > dfPixelSizeY)
        {
            nBadCount++;
        }
        nSamplePoints++;
    }

    return nBadCount == nSamplePoints;
}

/************************************************************************/
/*                      GDALSuggestedWarpOutput2()                      */
/************************************************************************/

/**
 * Suggest output file size.
 *
 * This function is used to suggest the size, and georeferenced extents
 * appropriate given the indicated transformation and input file.  It walks
 * the edges of the input file (approximately 20 sample points along each
 * edge) transforming into output coordinates in order to get an extents box.
 *
 * Then a resolution is computed with the intent that the length of the
 * distance from the top left corner of the output imagery to the bottom right
 * corner would represent the same number of pixels as in the source image.
 * Note that if the image is somewhat rotated the diagonal taken isn't of the
 * whole output bounding rectangle, but instead of the locations where the
 * top/left and bottom/right corners transform.  The output pixel size is
 * always square.  This is intended to approximately preserve the resolution
 * of the input data in the output file.
 *
 * The values returned in padfGeoTransformOut, pnPixels and pnLines are
 * the suggested number of pixels and lines for the output file, and the
 * geotransform relating those pixels to the output georeferenced coordinates.
 *
 * The trickiest part of using the function is ensuring that the
 * transformer created is from source file pixel/line coordinates to
 * output file georeferenced coordinates.  This can be accomplished with
 * GDALCreateGenImgProjTransformer() by passing a NULL for the hDstDS.
 *
 * @param hSrcDS the input image (it is assumed the whole input image is
 * being transformed).
 * @param pfnTransformer the transformer function.
 * @param pTransformArg the callback data for the transformer function.
 * @param padfGeoTransformOut the array of six doubles in which the suggested
 * geotransform is returned.
 * @param pnPixels int in which the suggest pixel width of output is returned.
 * @param pnLines int in which the suggest pixel height of output is returned.
 * @param padfExtent Four entry array to return extents as (xmin, ymin, xmax,
 * ymax).
 * @param nOptions Options flags. Zero or GDAL_SWO_ROUND_UP_SIZE  to ask *pnPixels
 * and *pnLines to be rounded up instead of being rounded to the closes integer, or
 * GDAL_SWO_FORCE_SQUARE_PIXEL to indicate that the generated pixel size is a square.
 *
 * @return CE_None if successful or CE_Failure otherwise.
 */

CPLErr CPL_STDCALL GDALSuggestedWarpOutput2(GDALDatasetH hSrcDS,
                                            GDALTransformerFunc pfnTransformer,
                                            void *pTransformArg,
                                            double *padfGeoTransformOut,
                                            int *pnPixels, int *pnLines,
                                            double *padfExtent, int nOptions)
{
    VALIDATE_POINTER1(hSrcDS, "GDALSuggestedWarpOutput2", CE_Failure);

    const bool bIsGDALGenImgProjTransform{
        pTransformArg &&
        GDALIsTransformer(pTransformArg, GDAL_GEN_IMG_TRANSFORMER_CLASS_NAME)};

    /* -------------------------------------------------------------------- */
    /*      Setup sample points all around the edge of the input raster.    */
    /* -------------------------------------------------------------------- */
    if (bIsGDALGenImgProjTransform)
    {
        // In case CHECK_WITH_INVERT_PROJ has been modified.
        GDALRefreshGenImgProjTransformer(pTransformArg);
    }
    else if (GDALIsTransformer(pTransformArg,
                               GDAL_APPROX_TRANSFORMER_CLASS_NAME))
    {
        // In case CHECK_WITH_INVERT_PROJ has been modified.
        GDALRefreshApproxTransformer(pTransformArg);
    }

    const int nInXSize = GDALGetRasterXSize(hSrcDS);
    const int nInYSize = GDALGetRasterYSize(hSrcDS);

    /* ------------------------------------------------------------- */
    /* Special case for warping on the same (or null) CRS.           */
    /* ------------------------------------------------------------- */
    if ((!nOptions || (nOptions & GDAL_SWO_FORCE_SQUARE_PIXEL) == 0) &&
        pTransformArg && bIsGDALGenImgProjTransform)
    {
        const GDALGenImgProjTransformInfo *psInfo =
            static_cast<const GDALGenImgProjTransformInfo *>(pTransformArg);

        if (!psInfo->sSrcParams.pTransformer &&
            !psInfo->bHasCustomTransformationPipeline &&
            !psInfo->sDstParams.pTransformer &&
            psInfo->sSrcParams.adfGeoTransform[2] == 0 &&
            psInfo->sSrcParams.adfGeoTransform[4] == 0 &&
            psInfo->sDstParams.adfGeoTransform[0] == 0 &&
            psInfo->sDstParams.adfGeoTransform[1] == 1 &&
            psInfo->sDstParams.adfGeoTransform[2] == 0 &&
            psInfo->sDstParams.adfGeoTransform[3] == 0 &&
            psInfo->sDstParams.adfGeoTransform[4] == 0 &&
            psInfo->sDstParams.adfGeoTransform[5] == 1)
        {
            const OGRSpatialReference *poSourceCRS = nullptr;
            const OGRSpatialReference *poTargetCRS = nullptr;

            if (psInfo->pReprojectArg)
            {
                const GDALReprojectionTransformInfo *psRTI =
                    static_cast<const GDALReprojectionTransformInfo *>(
                        psInfo->pReprojectArg);
                poSourceCRS = psRTI->poForwardTransform->GetSourceCS();
                poTargetCRS = psRTI->poForwardTransform->GetTargetCS();
            }

            if ((!poSourceCRS && !poTargetCRS) ||
                (poSourceCRS && poTargetCRS &&
                 poSourceCRS->IsSame(poTargetCRS)))
            {

                const bool bNorthUp{psInfo->sSrcParams.adfGeoTransform[5] <
                                    0.0};

                memcpy(padfGeoTransformOut, psInfo->sSrcParams.adfGeoTransform,
                       sizeof(double) * 6);

                if (!bNorthUp)
                {
                    padfGeoTransformOut[3] = padfGeoTransformOut[3] +
                                             nInYSize * padfGeoTransformOut[5];
                    padfGeoTransformOut[5] = -padfGeoTransformOut[5];
                }

                *pnPixels = nInXSize;
                *pnLines = nInYSize;

                // Calculate extent from hSrcDS
                if (padfExtent)
                {
                    padfExtent[0] = psInfo->sSrcParams.adfGeoTransform[0];
                    padfExtent[1] =
                        psInfo->sSrcParams.adfGeoTransform[3] +
                        nInYSize * psInfo->sSrcParams.adfGeoTransform[5];
                    padfExtent[2] =
                        psInfo->sSrcParams.adfGeoTransform[0] +
                        nInXSize * psInfo->sSrcParams.adfGeoTransform[1];
                    padfExtent[3] = psInfo->sSrcParams.adfGeoTransform[3];
                    if (!bNorthUp)
                    {
                        std::swap(padfExtent[1], padfExtent[3]);
                    }
                }
                return CE_None;
            }
        }
    }

    const int N_PIXELSTEP = 50;
    int nSteps = static_cast<int>(
        static_cast<double>(std::min(nInYSize, nInXSize)) / N_PIXELSTEP + 0.5);
    if (nSteps < 20)
        nSteps = 20;
    else if (nSteps > 100)
        nSteps = 100;

    // TODO(rouault): How is this goto retry supposed to work?  Added in r20537.
    // Does redoing the same malloc multiple times work?  If it is needed, can
    // it be converted to a tigher while loop around the MALLOC3s and free?  Is
    // the point to try with the full requested steps.  Then, if there is not
    // enough memory, back off and try with just 20 steps?
retry:
    int nStepsPlusOne = nSteps + 1;
    int nSampleMax = nStepsPlusOne * nStepsPlusOne;

    double dfStep = 1.0 / nSteps;
    double *padfY = nullptr;
    double *padfZ = nullptr;
    double *padfYRevert = nullptr;
    double *padfZRevert = nullptr;

    int *pabSuccess = static_cast<int *>(
        VSI_MALLOC3_VERBOSE(sizeof(int), nStepsPlusOne, nStepsPlusOne));
    double *padfX = static_cast<double *>(
        VSI_MALLOC3_VERBOSE(sizeof(double) * 3, nStepsPlusOne, nStepsPlusOne));
    double *padfXRevert = static_cast<double *>(
        VSI_MALLOC3_VERBOSE(sizeof(double) * 3, nStepsPlusOne, nStepsPlusOne));
    if (pabSuccess == nullptr || padfX == nullptr || padfXRevert == nullptr)
    {
        CPLFree(padfX);
        CPLFree(padfXRevert);
        CPLFree(pabSuccess);
        if (nSteps > 20)
        {
            nSteps = 20;
            goto retry;
        }
        return CE_Failure;
    }

    padfY = padfX + nSampleMax;
    padfZ = padfX + nSampleMax * 2;
    padfYRevert = padfXRevert + nSampleMax;
    padfZRevert = padfXRevert + nSampleMax * 2;

    // Take N_STEPS steps.
    for (int iStep = 0; iStep <= nSteps; iStep++)
    {
        double dfRatio = (iStep == nSteps) ? 1.0 : iStep * dfStep;
        int iStep2 = iStep;

        // Along top.
        padfX[iStep2] = dfRatio * nInXSize;
        padfY[iStep2] = 0.0;
        padfZ[iStep2] = 0.0;

        // Along bottom.
        iStep2 += nStepsPlusOne;
        padfX[iStep2] = dfRatio * nInXSize;
        padfY[iStep2] = nInYSize;
        padfZ[iStep2] = 0.0;

        // Along left.
        iStep2 += nStepsPlusOne;
        padfX[iStep2] = 0.0;
        padfY[iStep2] = dfRatio * nInYSize;
        padfZ[iStep2] = 0.0;

        // Along right.
        iStep2 += nStepsPlusOne;
        padfX[iStep2] = nInXSize;
        padfY[iStep2] = dfRatio * nInYSize;
        padfZ[iStep2] = 0.0;
    }

    int nSamplePoints = 4 * nStepsPlusOne;

    memset(pabSuccess, 1, sizeof(int) * nSampleMax);

    /* -------------------------------------------------------------------- */
    /*      Transform them to the output coordinate system.                 */
    /* -------------------------------------------------------------------- */
    {
        CPLTurnFailureIntoWarningBackuper oErrorsToWarnings{};
        pfnTransformer(pTransformArg, FALSE, nSamplePoints, padfX, padfY, padfZ,
                       pabSuccess);
    }
    constexpr int SIGN_FINAL_UNINIT = -2;
    constexpr int SIGN_FINAL_INVALID = 0;
    int iSignDiscontinuity = SIGN_FINAL_UNINIT;
    int nFailedCount = 0;
    const int iSignArray[2] = {-1, 1};
    for (int i = 0; i < nSamplePoints; i++)
    {
        if (pabSuccess[i])
        {
            // Fix for https://trac.osgeo.org/gdal/ticket/7243
            // where echo "-2050000.000 2050000.000" |
            //              gdaltransform -s_srs EPSG:3411 -t_srs EPSG:4326
            // gives "-180 63.691332898492"
            // but we would rather like 180
            if (iSignDiscontinuity == 1 || iSignDiscontinuity == -1)
            {
                if (!((iSignDiscontinuity * padfX[i] > 0 &&
                       iSignDiscontinuity * padfX[i] <= 180.0) ||
                      (fabs(padfX[i] - iSignDiscontinuity * -180.0) < 1e-8)))
                {
                    iSignDiscontinuity = SIGN_FINAL_INVALID;
                }
            }
            else if (iSignDiscontinuity == SIGN_FINAL_UNINIT)
            {
                for (const auto &iSign : iSignArray)
                {
                    if ((iSign * padfX[i] > 0 && iSign * padfX[i] <= 180.0) ||
                        (fabs(padfX[i] - iSign * -180.0) < 1e-8))
                    {
                        iSignDiscontinuity = iSign;
                        break;
                    }
                }
                if (iSignDiscontinuity == SIGN_FINAL_UNINIT)
                {
                    iSignDiscontinuity = SIGN_FINAL_INVALID;
                }
            }
        }
        else
        {
            nFailedCount++;
        }
    }

    if (iSignDiscontinuity == 1 || iSignDiscontinuity == -1)
    {
        for (int i = 0; i < nSamplePoints; i++)
        {
            if (pabSuccess[i])
            {
                if (fabs(padfX[i] - iSignDiscontinuity * -180.0) < 1e-8)
                {
                    double axTemp[2] = {iSignDiscontinuity * -180.0,
                                        iSignDiscontinuity * 180.0};
                    double ayTemp[2] = {padfY[i], padfY[i]};
                    double azTemp[2] = {padfZ[i], padfZ[i]};
                    int abSuccess[2] = {FALSE, FALSE};
                    CPLTurnFailureIntoWarningBackuper oErrorsToWarnings{};
                    if (pfnTransformer(pTransformArg, TRUE, 2, axTemp, ayTemp,
                                       azTemp, abSuccess) &&
                        fabs(axTemp[0] - axTemp[1]) < 1e-8 &&
                        fabs(ayTemp[0] - ayTemp[1]) < 1e-8)
                    {
                        padfX[i] = iSignDiscontinuity * 180.0;
                    }
                }
            }
        }
    }

    /* -------------------------------------------------------------------- */
    /*      Check if the computed target coordinates are revertable.        */
    /*      If not, try the detailed grid sampling.                         */
    /* -------------------------------------------------------------------- */
    if (nFailedCount)
    {
        CPLDebug("WARP", "At least one point failed after direct transform");
    }
    else
    {
        memcpy(padfXRevert, padfX, nSamplePoints * sizeof(double));
        memcpy(padfYRevert, padfY, nSamplePoints * sizeof(double));
        memcpy(padfZRevert, padfZ, nSamplePoints * sizeof(double));
        {
            CPLTurnFailureIntoWarningBackuper oErrorsToWarnings{};
            pfnTransformer(pTransformArg, TRUE, nSamplePoints, padfXRevert,
                           padfYRevert, padfZRevert, pabSuccess);
        }

        for (int i = 0; nFailedCount == 0 && i < nSamplePoints; i++)
        {
            if (!pabSuccess[i])
            {
                nFailedCount++;
                break;
            }

            double dfRatio = (i % nStepsPlusOne) * dfStep;
            if (dfRatio > 0.99)
                dfRatio = 1.0;

            double dfExpectedX = 0.0;
            double dfExpectedY = 0.0;
            if (i < nStepsPlusOne)
            {
                dfExpectedX = dfRatio * nInXSize;
            }
            else if (i < 2 * nStepsPlusOne)
            {
                dfExpectedX = dfRatio * nInXSize;
                dfExpectedY = nInYSize;
            }
            else if (i < 3 * nStepsPlusOne)
            {
                dfExpectedY = dfRatio * nInYSize;
            }
            else
            {
                dfExpectedX = nInXSize;
                dfExpectedY = dfRatio * nInYSize;
            }

            if (fabs(padfXRevert[i] - dfExpectedX) >
                    nInXSize / static_cast<double>(nSteps) ||
                fabs(padfYRevert[i] - dfExpectedY) >
                    nInYSize / static_cast<double>(nSteps))
                nFailedCount++;
        }
        if (nFailedCount != 0)
            CPLDebug("WARP",
                     "At least one point failed after revert transform");
    }

    /* -------------------------------------------------------------------- */
    /*      If any of the edge points failed to transform, we need to       */
    /*      build a fairly detailed internal grid of points instead to      */
    /*      help identify the area that is transformable.                   */
    /* -------------------------------------------------------------------- */
    if (nFailedCount)
    {
        nSamplePoints = 0;

        // Take N_STEPS steps.
        for (int iStep = 0; iStep <= nSteps; iStep++)
        {
            double dfRatio = (iStep == nSteps) ? 1.0 : iStep * dfStep;

            for (int iStep2 = 0; iStep2 <= nSteps; iStep2++)
            {
                const double dfRatio2 =
                    iStep2 == nSteps ? 1.0 : iStep2 * dfStep;

                // From top to bottom, from left to right.
                padfX[nSamplePoints] = dfRatio2 * nInXSize;
                padfY[nSamplePoints] = dfRatio * nInYSize;
                padfZ[nSamplePoints] = 0.0;
                nSamplePoints++;
            }
        }

        CPLAssert(nSamplePoints == nSampleMax);

        {
            CPLTurnFailureIntoWarningBackuper oErrorsToWarnings{};
            pfnTransformer(pTransformArg, FALSE, nSamplePoints, padfX, padfY,
                           padfZ, pabSuccess);
        }
    }

    /* -------------------------------------------------------------------- */
    /*      Collect the bounds, ignoring any failed points.                 */
    /* -------------------------------------------------------------------- */
    double dfMinXOut = 0.0;
    double dfMinYOut = 0.0;
    double dfMaxXOut = 0.0;
    double dfMaxYOut = 0.0;
    bool bGotInitialPoint = false;

    nFailedCount = 0;
    for (int i = 0; i < nSamplePoints; i++)
    {
        int x_i = 0;
        int y_i = 0;

        if (nSamplePoints == nSampleMax)
        {
            x_i = i % nStepsPlusOne;
            y_i = i / nStepsPlusOne;
        }
        else
        {
            if (i < 2 * nStepsPlusOne)
            {
                x_i = i % nStepsPlusOne;
                y_i = (i < nStepsPlusOne) ? 0 : nSteps;
            }
        }

        if (x_i > 0 && (pabSuccess[i - 1] || pabSuccess[i]))
        {
            double x_out_before = padfX[i - 1];
            double x_out_after = padfX[i];
            int nIter = 0;
            double x_in_before =
                static_cast<double>(x_i - 1) * nInXSize / nSteps;
            double x_in_after = static_cast<double>(x_i) * nInXSize / nSteps;
            int invalid_before = !(pabSuccess[i - 1]);
            int invalid_after = !(pabSuccess[i]);

            // Detect discontinuity in target coordinates when the target x
            // coordinates change sign. This may be a false positive when the
            // target tx is around 0 Dichotomic search to reduce the interval
            // to near the discontinuity and get a better out extent.
            while ((invalid_before || invalid_after ||
                    x_out_before * x_out_after < 0.0) &&
                   nIter < 16)
            {
                double x = (x_in_before + x_in_after) / 2.0;
                double y = static_cast<double>(y_i) * nInYSize / nSteps;
                double z = 0.0;
                int bSuccess = TRUE;
                if (pfnTransformer(pTransformArg, FALSE, 1, &x, &y, &z,
                                   &bSuccess) &&
                    bSuccess)
                {
                    if (bGotInitialPoint)
                    {
                        dfMinXOut = std::min(dfMinXOut, x);
                        dfMinYOut = std::min(dfMinYOut, y);
                        dfMaxXOut = std::max(dfMaxXOut, x);
                        dfMaxYOut = std::max(dfMaxYOut, y);
                    }
                    else
                    {
                        bGotInitialPoint = true;
                        dfMinXOut = x;
                        dfMaxXOut = x;
                        dfMinYOut = y;
                        dfMaxYOut = y;
                    }

                    if (invalid_before || x_out_before * x < 0)
                    {
                        invalid_after = FALSE;
                        x_in_after = (x_in_before + x_in_after) / 2.0;
                        x_out_after = x;
                    }
                    else
                    {
                        invalid_before = FALSE;
                        x_out_before = x;
                        x_in_before = (x_in_before + x_in_after) / 2.0;
                    }
                }
                else
                {
                    if (invalid_before)
                    {
                        x_in_before = (x_in_before + x_in_after) / 2.0;
                    }
                    else if (invalid_after)
                    {
                        x_in_after = (x_in_before + x_in_after) / 2.0;
                    }
                    else
                    {
                        break;
                    }
                }
                nIter++;
            }
        }

        if (!pabSuccess[i])
        {
            nFailedCount++;
            continue;
        }

        if (bGotInitialPoint)
        {
            dfMinXOut = std::min(dfMinXOut, padfX[i]);
            dfMinYOut = std::min(dfMinYOut, padfY[i]);
            dfMaxXOut = std::max(dfMaxXOut, padfX[i]);
            dfMaxYOut = std::max(dfMaxYOut, padfY[i]);
        }
        else
        {
            bGotInitialPoint = true;
            dfMinXOut = padfX[i];
            dfMaxXOut = padfX[i];
            dfMinYOut = padfY[i];
            dfMaxYOut = padfY[i];
        }
    }

    if (nFailedCount > nSamplePoints - 10)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Too many points (%d out of %d) failed to transform, "
                 "unable to compute output bounds.",
                 nFailedCount, nSamplePoints);

        CPLFree(padfX);
        CPLFree(padfXRevert);
        CPLFree(pabSuccess);

        return CE_Failure;
    }

    if (nFailedCount)
        CPLDebug("GDAL",
                 "GDALSuggestedWarpOutput(): %d out of %d points failed to "
                 "transform.",
                 nFailedCount, nSamplePoints);

    bool bIsGeographicCoordsDeg = false;
    if (bIsGDALGenImgProjTransform)
    {
        const GDALGenImgProjTransformInfo *pGIPTI =
            static_cast<const GDALGenImgProjTransformInfo *>(pTransformArg);
        if (pGIPTI->sSrcParams.pTransformer == GDALGeoLocTransform &&
            pGIPTI->sDstParams.pTransformer == nullptr &&
            pGIPTI->sDstParams.adfGeoTransform[0] == 0 &&
            pGIPTI->sDstParams.adfGeoTransform[1] == 1 &&
            pGIPTI->sDstParams.adfGeoTransform[2] == 0 &&
            pGIPTI->sDstParams.adfGeoTransform[3] == 0 &&
            pGIPTI->sDstParams.adfGeoTransform[4] == 0 &&
            pGIPTI->sDstParams.adfGeoTransform[5] == 1)
        {
            /* --------------------------------------------------------------------
             */
            /*      Special case for geolocation array, to quickly find the
             * bounds. */
            /* --------------------------------------------------------------------
             */
            const GDALGeoLocTransformInfo *pGLTI =
                static_cast<const GDALGeoLocTransformInfo *>(
                    pGIPTI->sSrcParams.pTransformArg);

            if (pGIPTI->pReproject == nullptr)
            {
                const char *pszGLSRS =
                    CSLFetchNameValue(pGLTI->papszGeolocationInfo, "SRS");
                if (pszGLSRS == nullptr)
                {
                    bIsGeographicCoordsDeg = true;
                }
                else
                {
                    OGRSpatialReference oSRS;
                    if (oSRS.SetFromUserInput(pszGLSRS) == OGRERR_NONE &&
                        oSRS.IsGeographic() &&
                        std::fabs(oSRS.GetAngularUnits() -
                                  CPLAtof(SRS_UA_DEGREE_CONV)) < 1e-9)
                    {
                        bIsGeographicCoordsDeg = true;
                    }
                }
            }

            for (const auto &xy :
                 {std::pair<double, double>(pGLTI->dfMinX, pGLTI->dfYAtMinX),
                  std::pair<double, double>(pGLTI->dfXAtMinY, pGLTI->dfMinY),
                  std::pair<double, double>(pGLTI->dfMaxX, pGLTI->dfYAtMaxX),
                  std::pair<double, double>(pGLTI->dfXAtMaxY, pGLTI->dfMaxY)})
            {
                double x = xy.first;
                double y = xy.second;
                if (pGLTI->bSwapXY)
                {
                    std::swap(x, y);
                }
                double xOut = std::numeric_limits<double>::quiet_NaN();
                double yOut = std::numeric_limits<double>::quiet_NaN();
                if (pGIPTI->pReproject == nullptr ||
                    pGIPTI->pReproject(pGIPTI->pReprojectArg, false, 1, &x, &y,
                                       nullptr, nullptr))
                {
                    xOut = x;
                    yOut = y;
                }
                dfMinXOut = std::min(dfMinXOut, xOut);
                dfMinYOut = std::min(dfMinYOut, yOut);
                dfMaxXOut = std::max(dfMaxXOut, xOut);
                dfMaxYOut = std::max(dfMaxYOut, yOut);
            }
        }
        else if (pGIPTI->sSrcParams.pTransformer == nullptr &&
                 pGIPTI->sDstParams.pTransformer == nullptr &&
                 pGIPTI->pReproject == GDALReprojectionTransform &&
                 pGIPTI->sDstParams.adfGeoTransform[0] == 0 &&
                 pGIPTI->sDstParams.adfGeoTransform[1] == 1 &&
                 pGIPTI->sDstParams.adfGeoTransform[2] == 0 &&
                 pGIPTI->sDstParams.adfGeoTransform[3] == 0 &&
                 pGIPTI->sDstParams.adfGeoTransform[4] == 0 &&
                 pGIPTI->sDstParams.adfGeoTransform[5] == 1)
        {
            /* ------------------------------------------------------------- */
            /* Special case for warping using source geotransform and        */
            /* reprojection to deal with the poles.                          */
            /* ------------------------------------------------------------- */
            const GDALReprojectionTransformInfo *psRTI =
                static_cast<const GDALReprojectionTransformInfo *>(
                    pGIPTI->pReprojectArg);
            const OGRSpatialReference *poSourceCRS =
                psRTI->poForwardTransform->GetSourceCS();
            const OGRSpatialReference *poTargetCRS =
                psRTI->poForwardTransform->GetTargetCS();
            if (poTargetCRS != nullptr &&
                psRTI->poReverseTransform != nullptr &&
                poTargetCRS->IsGeographic() &&
                fabs(poTargetCRS->GetAngularUnits() -
                     CPLAtof(SRS_UA_DEGREE_CONV)) < 1e-9 &&
                (!poSourceCRS || !poSourceCRS->IsGeographic()))
            {
                bIsGeographicCoordsDeg = true;

                std::unique_ptr<CPLConfigOptionSetter> poSetter;
                if (pGIPTI->bCheckWithInvertPROJ)
                {
                    // CHECK_WITH_INVERT_PROJ=YES prevent reliable
                    // transformation of poles.
                    poSetter = std::make_unique<CPLConfigOptionSetter>(
                        "CHECK_WITH_INVERT_PROJ", "NO", false);
                    GDALRefreshGenImgProjTransformer(pTransformArg);
                    // GDALRefreshGenImgProjTransformer() has invalidated psRTI
                    psRTI = static_cast<const GDALReprojectionTransformInfo *>(
                        pGIPTI->pReprojectArg);
                }

                for (const auto &sign : iSignArray)
                {
                    double X = 0.0;
                    const double Yinit = 90.0 * sign;
                    double Y = Yinit;
                    if (psRTI->poReverseTransform->Transform(1, &X, &Y))
                    {
                        const auto invGT =
                            pGIPTI->sSrcParams.adfInvGeoTransform;
                        const double x = invGT[0] + X * invGT[1] + Y * invGT[2];
                        const double y = invGT[3] + X * invGT[4] + Y * invGT[5];
                        constexpr double EPSILON = 1e-5;
                        if (x >= -EPSILON && x <= nInXSize + EPSILON &&
                            y >= -EPSILON && y <= nInYSize + EPSILON)
                        {
                            if (psRTI->poForwardTransform->Transform(1, &X,
                                                                     &Y) &&
                                fabs(Y - Yinit) <= 1e-6)
                            {
                                bool bMinXMaxXSet = false;
                                if (poSourceCRS)
                                {
                                    const char *pszProjection =
                                        poSourceCRS->GetAttrValue("PROJECTION");
                                    if (pszProjection &&
                                        EQUAL(pszProjection,
                                              SRS_PT_ORTHOGRAPHIC))
                                    {
                                        const double dfLon0 =
                                            poSourceCRS->GetNormProjParm(
                                                SRS_PP_CENTRAL_MERIDIAN, 0.0);
                                        dfMinXOut = dfLon0 - 90;
                                        dfMaxXOut = dfLon0 + 90;
                                        bMinXMaxXSet = true;
                                    }
                                }
                                if (!bMinXMaxXSet)
                                {
                                    dfMinXOut = -180;
                                    dfMaxXOut = 180;
                                }
                                if (sign < 0)
                                    dfMinYOut = Yinit;
                                else
                                    dfMaxYOut = Yinit;
                            }
                        }
                    }
                }

                if (poSetter)
                {
                    poSetter.reset();
                    GDALRefreshGenImgProjTransformer(pTransformArg);
                    pGIPTI = static_cast<const GDALGenImgProjTransformInfo *>(
                        pTransformArg);
                    psRTI = static_cast<const GDALReprojectionTransformInfo *>(
                        pGIPTI->pReprojectArg);
                    poSourceCRS = psRTI->poForwardTransform->GetSourceCS();
                    poTargetCRS = psRTI->poForwardTransform->GetTargetCS();
                }
            }

            // Use TransformBounds() to handle more particular cases
            if (poSourceCRS != nullptr && poTargetCRS != nullptr &&
                pGIPTI->sSrcParams.adfGeoTransform[1] != 0 &&
                pGIPTI->sSrcParams.adfGeoTransform[2] == 0 &&
                pGIPTI->sSrcParams.adfGeoTransform[4] == 0 &&
                pGIPTI->sSrcParams.adfGeoTransform[5] != 0)
            {
                const double dfULX = pGIPTI->sSrcParams.adfGeoTransform[0];
                const double dfULY = pGIPTI->sSrcParams.adfGeoTransform[3];
                const double dfLRX =
                    dfULX + pGIPTI->sSrcParams.adfGeoTransform[1] * nInXSize;
                const double dfLRY =
                    dfULY + pGIPTI->sSrcParams.adfGeoTransform[5] * nInYSize;
                const double dfMinSrcX = std::min(dfULX, dfLRX);
                const double dfMinSrcY = std::min(dfULY, dfLRY);
                const double dfMaxSrcX = std::max(dfULX, dfLRX);
                const double dfMaxSrcY = std::max(dfULY, dfLRY);
                double dfTmpMinXOut = std::numeric_limits<double>::max();
                double dfTmpMinYOut = std::numeric_limits<double>::max();
                double dfTmpMaxXOut = std::numeric_limits<double>::min();
                double dfTmpMaxYOut = std::numeric_limits<double>::min();
                if (psRTI->poForwardTransform->TransformBounds(
                        dfMinSrcX, dfMinSrcY, dfMaxSrcX, dfMaxSrcY,
                        &dfTmpMinXOut, &dfTmpMinYOut, &dfTmpMaxXOut,
                        &dfTmpMaxYOut,
                        2))  // minimum number of points as we already have a
                             // logic above to sample
                {
                    dfMinXOut = std::min(dfMinXOut, dfTmpMinXOut);
                    dfMinYOut = std::min(dfMinYOut, dfTmpMinYOut);
                    dfMaxXOut = std::max(dfMaxXOut, dfTmpMaxXOut);
                    dfMaxYOut = std::max(dfMaxYOut, dfTmpMaxYOut);
                }
            }
        }
    }

    /* -------------------------------------------------------------------- */
    /*      Compute the distance in "georeferenced" units from the top      */
    /*      corner of the transformed input image to the bottom left        */
    /*      corner of the transformed input.  Use this distance to          */
    /*      compute an approximate pixel size in the output                 */
    /*      georeferenced coordinates.                                      */
    /* -------------------------------------------------------------------- */
    double dfDiagonalDist = 0.0;
    double dfDeltaX = 0.0;
    double dfDeltaY = 0.0;

    if (pabSuccess[0] && pabSuccess[nSamplePoints - 1])
    {
        dfDeltaX = padfX[nSamplePoints - 1] - padfX[0];
        dfDeltaY = padfY[nSamplePoints - 1] - padfY[0];
        // In some cases this can result in 0 values. See #5980
        // Fallback to safer method in that case.
    }
    if (dfDeltaX == 0.0 || dfDeltaY == 0.0)
    {
        dfDeltaX = dfMaxXOut - dfMinXOut;
        dfDeltaY = dfMaxYOut - dfMinYOut;
    }

    dfDiagonalDist = sqrt(dfDeltaX * dfDeltaX + dfDeltaY * dfDeltaY);

    /* -------------------------------------------------------------------- */
    /*      Compute a pixel size from this.                                 */
    /* -------------------------------------------------------------------- */
    const double dfPixelSize =
        dfDiagonalDist / sqrt(static_cast<double>(nInXSize) * nInXSize +
                              static_cast<double>(nInYSize) * nInYSize);

    const double dfPixels = (dfMaxXOut - dfMinXOut) / dfPixelSize;
    const double dfLines = (dfMaxYOut - dfMinYOut) / dfPixelSize;

    const int knIntMaxMinusOne = std::numeric_limits<int>::max() - 1;
    if (dfPixels > knIntMaxMinusOne || dfLines > knIntMaxMinusOne)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Computed dimensions are too big : %.0f x %.0f",
                 dfPixels + 0.5, dfLines + 0.5);

        CPLFree(padfX);
        CPLFree(padfXRevert);
        CPLFree(pabSuccess);

        return CE_Failure;
    }

    if ((nOptions & GDAL_SWO_ROUND_UP_SIZE) != 0)
    {
        constexpr double EPS = 1e-5;
        *pnPixels = static_cast<int>(std::ceil(dfPixels - EPS));
        *pnLines = static_cast<int>(std::ceil(dfLines - EPS));
    }
    else
    {
        *pnPixels = static_cast<int>(dfPixels + 0.5);
        *pnLines = static_cast<int>(dfLines + 0.5);
    }

    double dfPixelSizeX = dfPixelSize;
    double dfPixelSizeY = dfPixelSize;

    const double adfRatioArray[] = {0.000, 0.001, 0.010, 0.100, 1.000};

    /* -------------------------------------------------------------------- */
    /*      Check that the right border is not completely out of source     */
    /*      image. If so, adjust the x pixel size a bit in the hope it will */
    /*      fit.                                                            */
    /* -------------------------------------------------------------------- */
    for (const auto &dfRatio : adfRatioArray)
    {
        const double dfTryPixelSizeX =
            dfPixelSizeX - dfPixelSizeX * dfRatio / *pnPixels;
        double adfExtent[4] = {dfMinXOut, dfMaxYOut - (*pnLines) * dfPixelSizeY,
                               dfMinXOut + (*pnPixels) * dfTryPixelSizeX,
                               dfMaxYOut};
        if (!GDALSuggestedWarpOutput2_MustAdjustForRightBorder(
                pfnTransformer, pTransformArg, adfExtent, *pnPixels, *pnLines,
                dfTryPixelSizeX, dfPixelSizeY))
        {
            dfPixelSizeX = dfTryPixelSizeX;
            break;
        }
    }

    /* -------------------------------------------------------------------- */
    /*      Check that the bottom border is not completely out of source    */
    /*      image. If so, adjust the y pixel size a bit in the hope it will */
    /*      fit.                                                            */
    /* -------------------------------------------------------------------- */
    for (const auto &dfRatio : adfRatioArray)
    {
        const double dfTryPixelSizeY =
            dfPixelSizeY - dfPixelSizeY * dfRatio / *pnLines;
        double adfExtent[4] = {
            dfMinXOut, dfMaxYOut - (*pnLines) * dfTryPixelSizeY,
            dfMinXOut + (*pnPixels) * dfPixelSizeX, dfMaxYOut};
        if (!GDALSuggestedWarpOutput2_MustAdjustForBottomBorder(
                pfnTransformer, pTransformArg, adfExtent, *pnPixels, *pnLines,
                dfPixelSizeX, dfTryPixelSizeY))
        {
            dfPixelSizeY = dfTryPixelSizeY;
            break;
        }
    }

    /* -------------------------------------------------------------------- */
    /*      Recompute some bounds so that all return values are consistent  */
    /* -------------------------------------------------------------------- */
    double dfMaxXOutNew = dfMinXOut + (*pnPixels) * dfPixelSizeX;
    if (bIsGeographicCoordsDeg &&
        ((dfMaxXOut <= 180 && dfMaxXOutNew > 180) || dfMaxXOut == 180))
    {
        dfMaxXOut = 180;
        dfPixelSizeX = (dfMaxXOut - dfMinXOut) / *pnPixels;
    }
    else
    {
        dfMaxXOut = dfMaxXOutNew;
    }

    double dfMinYOutNew = dfMaxYOut - (*pnLines) * dfPixelSizeY;
    if (bIsGeographicCoordsDeg && dfMinYOut >= -90 && dfMinYOutNew < -90)
    {
        dfMinYOut = -90;
        dfPixelSizeY = (dfMaxYOut - dfMinYOut) / *pnLines;
    }
    else
    {
        dfMinYOut = dfMinYOutNew;
    }

    /* -------------------------------------------------------------------- */
    /*      Return raw extents.                                             */
    /* -------------------------------------------------------------------- */
    padfExtent[0] = dfMinXOut;
    padfExtent[1] = dfMinYOut;
    padfExtent[2] = dfMaxXOut;
    padfExtent[3] = dfMaxYOut;

    /* -------------------------------------------------------------------- */
    /*      Set the output geotransform.                                    */
    /* -------------------------------------------------------------------- */
    padfGeoTransformOut[0] = dfMinXOut;
    padfGeoTransformOut[1] = dfPixelSizeX;
    padfGeoTransformOut[2] = 0.0;
    padfGeoTransformOut[3] = dfMaxYOut;
    padfGeoTransformOut[4] = 0.0;
    padfGeoTransformOut[5] = -dfPixelSizeY;

    CPLFree(padfX);
    CPLFree(padfXRevert);
    CPLFree(pabSuccess);

    return CE_None;
}

/************************************************************************/
/*                    GetCurrentCheckWithInvertPROJ()                   */
/************************************************************************/

static bool GetCurrentCheckWithInvertPROJ()
{
    return CPLTestBool(CPLGetConfigOption("CHECK_WITH_INVERT_PROJ", "NO"));
}

/************************************************************************/
/*               GDALCreateGenImgProjTransformerInternal()              */
/************************************************************************/

static void *GDALCreateSimilarGenImgProjTransformer(void *hTransformArg,
                                                    double dfRatioX,
                                                    double dfRatioY);

static GDALGenImgProjTransformInfo *GDALCreateGenImgProjTransformerInternal()
{
    /* -------------------------------------------------------------------- */
    /*      Initialize the transform info.                                  */
    /* -------------------------------------------------------------------- */
    GDALGenImgProjTransformInfo *psInfo =
        static_cast<GDALGenImgProjTransformInfo *>(
            CPLCalloc(sizeof(GDALGenImgProjTransformInfo), 1));

    memcpy(psInfo->sTI.abySignature, GDAL_GTI2_SIGNATURE,
           strlen(GDAL_GTI2_SIGNATURE));
    psInfo->sTI.pszClassName = GDAL_GEN_IMG_TRANSFORMER_CLASS_NAME;
    psInfo->sTI.pfnTransform = GDALGenImgProjTransform;
    psInfo->sTI.pfnCleanup = GDALDestroyGenImgProjTransformer;
    psInfo->sTI.pfnSerialize = GDALSerializeGenImgProjTransformer;
    psInfo->sTI.pfnCreateSimilar = GDALCreateSimilarGenImgProjTransformer;

    psInfo->bCheckWithInvertPROJ = GetCurrentCheckWithInvertPROJ();
    psInfo->bHasCustomTransformationPipeline = false;

    return psInfo;
}

/************************************************************************/
/*                GDALCreateSimilarGenImgProjTransformer()              */
/************************************************************************/

static void *GDALCreateSimilarGenImgProjTransformer(void *hTransformArg,
                                                    double dfRatioX,
                                                    double dfRatioY)
{
    VALIDATE_POINTER1(hTransformArg, "GDALCreateSimilarGenImgProjTransformer",
                      nullptr);

    GDALGenImgProjTransformInfo *psInfo =
        static_cast<GDALGenImgProjTransformInfo *>(hTransformArg);

    GDALGenImgProjTransformInfo *psClonedInfo =
        GDALCreateGenImgProjTransformerInternal();

    memcpy(psClonedInfo, psInfo, sizeof(GDALGenImgProjTransformInfo));

    psClonedInfo->bCheckWithInvertPROJ = GetCurrentCheckWithInvertPROJ();

    if (psClonedInfo->sSrcParams.pTransformArg)
        psClonedInfo->sSrcParams.pTransformArg = GDALCreateSimilarTransformer(
            psInfo->sSrcParams.pTransformArg, dfRatioX, dfRatioY);
    else if (dfRatioX != 1.0 || dfRatioY != 1.0)
    {
        if (psClonedInfo->sSrcParams.adfGeoTransform[2] == 0.0 &&
            psClonedInfo->sSrcParams.adfGeoTransform[4] == 0.0)
        {
            psClonedInfo->sSrcParams.adfGeoTransform[1] *= dfRatioX;
            psClonedInfo->sSrcParams.adfGeoTransform[5] *= dfRatioY;
        }
        else
        {
            // If the x and y ratios are not equal, then we cannot really
            // compute a geotransform.
            psClonedInfo->sSrcParams.adfGeoTransform[1] *= dfRatioX;
            psClonedInfo->sSrcParams.adfGeoTransform[2] *= dfRatioX;
            psClonedInfo->sSrcParams.adfGeoTransform[4] *= dfRatioX;
            psClonedInfo->sSrcParams.adfGeoTransform[5] *= dfRatioX;
        }
        if (!GDALInvGeoTransform(psClonedInfo->sSrcParams.adfGeoTransform,
                                 psClonedInfo->sSrcParams.adfInvGeoTransform))
        {
            CPLError(CE_Failure, CPLE_AppDefined, "Cannot invert geotransform");
            GDALDestroyGenImgProjTransformer(psClonedInfo);
            return nullptr;
        }
    }

    if (psClonedInfo->pReprojectArg)
        psClonedInfo->pReprojectArg =
            GDALCloneTransformer(psInfo->pReprojectArg);

    if (psClonedInfo->sDstParams.pTransformArg)
        psClonedInfo->sDstParams.pTransformArg =
            GDALCloneTransformer(psInfo->sDstParams.pTransformArg);

    return psClonedInfo;
}

/************************************************************************/
/*                  GDALCreateGenImgProjTransformer()                   */
/************************************************************************/

/**
 * Create image to image transformer.
 *
 * This function creates a transformation object that maps from pixel/line
 * coordinates on one image to pixel/line coordinates on another image.  The
 * images may potentially be georeferenced in different coordinate systems,
 * and may used GCPs to map between their pixel/line coordinates and
 * georeferenced coordinates (as opposed to the default assumption that their
 * geotransform should be used).
 *
 * This transformer potentially performs three concatenated transformations.
 *
 * The first stage is from source image pixel/line coordinates to source
 * image georeferenced coordinates, and may be done using the geotransform,
 * or if not defined using a polynomial model derived from GCPs.  If GCPs
 * are used this stage is accomplished using GDALGCPTransform().
 *
 * The second stage is to change projections from the source coordinate system
 * to the destination coordinate system, assuming they differ.  This is
 * accomplished internally using GDALReprojectionTransform().
 *
 * The third stage is converting from destination image georeferenced
 * coordinates to destination image coordinates.  This is done using the
 * destination image geotransform, or if not available, using a polynomial
 * model derived from GCPs. If GCPs are used this stage is accomplished using
 * GDALGCPTransform().  This stage is skipped if hDstDS is NULL when the
 * transformation is created.
 *
 * @param hSrcDS source dataset, or NULL.
 * @param pszSrcWKT the coordinate system for the source dataset.  If NULL,
 * it will be read from the dataset itself.
 * @param hDstDS destination dataset (or NULL).
 * @param pszDstWKT the coordinate system for the destination dataset.  If
 * NULL, and hDstDS not NULL, it will be read from the destination dataset.
 * @param bGCPUseOK TRUE if GCPs should be used if the geotransform is not
 * available on the source dataset (not destination).
 * @param dfGCPErrorThreshold ignored/deprecated.
 * @param nOrder the maximum order to use for GCP derived polynomials if
 * possible.  Use 0 to autoselect, or -1 for thin plate splines.
 *
 * @return handle suitable for use GDALGenImgProjTransform(), and to be
 * deallocated with GDALDestroyGenImgProjTransformer().
 */

void *GDALCreateGenImgProjTransformer(GDALDatasetH hSrcDS,
                                      const char *pszSrcWKT,
                                      GDALDatasetH hDstDS,
                                      const char *pszDstWKT, int bGCPUseOK,
                                      CPL_UNUSED double dfGCPErrorThreshold,
                                      int nOrder)
{
    char **papszOptions = nullptr;

    if (pszSrcWKT != nullptr)
        papszOptions = CSLSetNameValue(papszOptions, "SRC_SRS", pszSrcWKT);
    if (pszDstWKT != nullptr)
        papszOptions = CSLSetNameValue(papszOptions, "DST_SRS", pszDstWKT);
    if (!bGCPUseOK)
        papszOptions = CSLSetNameValue(papszOptions, "GCPS_OK", "FALSE");
    if (nOrder != 0)
        papszOptions = CSLSetNameValue(papszOptions, "MAX_GCP_ORDER",
                                       CPLString().Printf("%d", nOrder));

    void *pRet = GDALCreateGenImgProjTransformer2(hSrcDS, hDstDS, papszOptions);
    CSLDestroy(papszOptions);

    return pRet;
}

/************************************************************************/
/*                          InsertCenterLong()                          */
/*                                                                      */
/*      Insert a CENTER_LONG Extension entry on a GEOGCS to indicate    */
/*      the center longitude of the dataset for wrapping purposes.      */
/************************************************************************/

static void InsertCenterLong(GDALDatasetH hDS, OGRSpatialReference *poSRS,
                             CPLStringList &aosOptions)

{
    if (!poSRS->IsGeographic() || std::fabs(poSRS->GetAngularUnits() -
                                            CPLAtof(SRS_UA_DEGREE_CONV)) > 1e-9)
    {
        return;
    }

    if (poSRS->GetExtension(nullptr, "CENTER_LONG"))
        return;

    /* -------------------------------------------------------------------- */
    /*      For now we only do this if we have a geotransform since         */
    /*      other forms require a bunch of extra work.                      */
    /* -------------------------------------------------------------------- */
    double adfGeoTransform[6] = {};

    if (GDALGetGeoTransform(hDS, adfGeoTransform) != CE_None)
        return;

    /* -------------------------------------------------------------------- */
    /*      Compute min/max longitude based on testing the four corners.    */
    /* -------------------------------------------------------------------- */
    const int nXSize = GDALGetRasterXSize(hDS);
    const int nYSize = GDALGetRasterYSize(hDS);

    const double dfMinLong =
        std::min(std::min(adfGeoTransform[0] + 0 * adfGeoTransform[1] +
                              0 * adfGeoTransform[2],
                          adfGeoTransform[0] + nXSize * adfGeoTransform[1] +
                              0 * adfGeoTransform[2]),
                 std::min(adfGeoTransform[0] + 0 * adfGeoTransform[1] +
                              nYSize * adfGeoTransform[2],
                          adfGeoTransform[0] + nXSize * adfGeoTransform[1] +
                              nYSize * adfGeoTransform[2]));
    const double dfMaxLong =
        std::max(std::max(adfGeoTransform[0] + 0 * adfGeoTransform[1] +
                              0 * adfGeoTransform[2],
                          adfGeoTransform[0] + nXSize * adfGeoTransform[1] +
                              0 * adfGeoTransform[2]),
                 std::max(adfGeoTransform[0] + 0 * adfGeoTransform[1] +
                              nYSize * adfGeoTransform[2],
                          adfGeoTransform[0] + nXSize * adfGeoTransform[1] +
                              nYSize * adfGeoTransform[2]));

    const double dfEpsilon =
        std::max(std::fabs(adfGeoTransform[1]), std::fabs(adfGeoTransform[2]));
    // If the raster covers more than 360 degree (allow an extra pixel),
    // give up
    constexpr double RELATIVE_EPSILON = 0.05;  // for numeric precision issues
    if (dfMaxLong - dfMinLong > 360.0 + dfEpsilon * (1 + RELATIVE_EPSILON))
        return;

    /* -------------------------------------------------------------------- */
    /*      Insert center long.                                             */
    /* -------------------------------------------------------------------- */
    const double dfCenterLong = (dfMaxLong + dfMinLong) / 2.0;
    aosOptions.SetNameValue("CENTER_LONG", CPLSPrintf("%g", dfCenterLong));
}

/************************************************************************/
/*                      GDALComputeAreaOfInterest()                     */
/************************************************************************/

bool GDALComputeAreaOfInterest(const OGRSpatialReference *poSRS,
                               double adfGT[6], int nXSize, int nYSize,
                               double &dfWestLongitudeDeg,
                               double &dfSouthLatitudeDeg,
                               double &dfEastLongitudeDeg,
                               double &dfNorthLatitudeDeg)
{
    bool ret = false;

    if (!poSRS)
        return false;

    OGRSpatialReference oSrcSRSHoriz(*poSRS);
    if (oSrcSRSHoriz.IsCompound())
    {
        oSrcSRSHoriz.StripVertical();
    }

    OGRSpatialReference *poGeog = oSrcSRSHoriz.CloneGeogCS();
    if (poGeog)
    {
        poGeog->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        poGeog->SetAngularUnits(SRS_UA_DEGREE, CPLAtof(SRS_UA_DEGREE_CONV));

        auto poCT = OGRCreateCoordinateTransformation(&oSrcSRSHoriz, poGeog);
        if (poCT)
        {
            poCT->SetEmitErrors(false);

            double x[4], y[4];
            x[0] = adfGT[0];
            y[0] = adfGT[3];
            x[1] = adfGT[0] + nXSize * adfGT[1];
            y[1] = adfGT[3];
            x[2] = adfGT[0];
            y[2] = adfGT[3] + nYSize * adfGT[5];
            x[3] = x[1];
            y[3] = y[2];
            int validity[4] = {false, false, false, false};
            poCT->Transform(4, x, y, nullptr, validity);
            dfWestLongitudeDeg = std::numeric_limits<double>::max();
            dfSouthLatitudeDeg = std::numeric_limits<double>::max();
            dfEastLongitudeDeg = -std::numeric_limits<double>::max();
            dfNorthLatitudeDeg = -std::numeric_limits<double>::max();
            for (int i = 0; i < 4; i++)
            {
                if (validity[i])
                {
                    ret = true;
                    dfWestLongitudeDeg = std::min(dfWestLongitudeDeg, x[i]);
                    dfSouthLatitudeDeg = std::min(dfSouthLatitudeDeg, y[i]);
                    dfEastLongitudeDeg = std::max(dfEastLongitudeDeg, x[i]);
                    dfNorthLatitudeDeg = std::max(dfNorthLatitudeDeg, y[i]);
                }
            }
            if (validity[0] && validity[1] && x[0] > x[1])
            {
                dfWestLongitudeDeg = x[0];
                dfEastLongitudeDeg = x[1];
            }
            if (ret && std::fabs(dfWestLongitudeDeg) <= 180 &&
                std::fabs(dfEastLongitudeDeg) <= 180 &&
                std::fabs(dfSouthLatitudeDeg) <= 90 &&
                std::fabs(dfNorthLatitudeDeg) <= 90)
            {
                CPLDebug("GDAL", "Computing area of interest: %g, %g, %g, %g",
                         dfWestLongitudeDeg, dfSouthLatitudeDeg,
                         dfEastLongitudeDeg, dfNorthLatitudeDeg);
            }
            else
            {
                CPLDebug("GDAL", "Could not compute area of interest");
                dfWestLongitudeDeg = 0;
                dfSouthLatitudeDeg = 0;
                dfEastLongitudeDeg = 0;
                dfNorthLatitudeDeg = 0;
            }
            OGRCoordinateTransformation::DestroyCT(poCT);
        }

        delete poGeog;
    }

    return ret;
}

bool GDALComputeAreaOfInterest(const OGRSpatialReference *poSRS, double dfX1,
                               double dfY1, double dfX2, double dfY2,
                               double &dfWestLongitudeDeg,
                               double &dfSouthLatitudeDeg,
                               double &dfEastLongitudeDeg,
                               double &dfNorthLatitudeDeg)
{
    bool ret = false;

    if (!poSRS)
        return false;

    OGRSpatialReference oSrcSRSHoriz(*poSRS);
    if (oSrcSRSHoriz.IsCompound())
    {
        oSrcSRSHoriz.StripVertical();
    }

    OGRSpatialReference *poGeog = oSrcSRSHoriz.CloneGeogCS();
    if (poGeog)
    {
        poGeog->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        auto poCT = OGRCreateCoordinateTransformation(&oSrcSRSHoriz, poGeog);
        if (poCT)
        {
            double x[4], y[4];
            x[0] = dfX1;
            y[0] = dfY1;
            x[1] = dfX2;
            y[1] = dfY1;
            x[2] = dfX1;
            y[2] = dfY2;
            x[3] = dfX2;
            y[3] = dfY2;
            int validity[4] = {false, false, false, false};
            poCT->Transform(4, x, y, nullptr, validity);
            dfWestLongitudeDeg = std::numeric_limits<double>::max();
            dfSouthLatitudeDeg = std::numeric_limits<double>::max();
            dfEastLongitudeDeg = -std::numeric_limits<double>::max();
            dfNorthLatitudeDeg = -std::numeric_limits<double>::max();
            for (int i = 0; i < 4; i++)
            {
                if (validity[i])
                {
                    ret = true;
                    dfWestLongitudeDeg = std::min(dfWestLongitudeDeg, x[i]);
                    dfSouthLatitudeDeg = std::min(dfSouthLatitudeDeg, y[i]);
                    dfEastLongitudeDeg = std::max(dfEastLongitudeDeg, x[i]);
                    dfNorthLatitudeDeg = std::max(dfNorthLatitudeDeg, y[i]);
                }
            }
            if (validity[0] && validity[1] && (dfX1 - dfX2) * (x[0] - x[1]) < 0)
            {
                dfWestLongitudeDeg = x[0];
                dfEastLongitudeDeg = x[1];
            }
            if (ret)
            {
                CPLDebug("GDAL", "Computing area of interest: %g, %g, %g, %g",
                         dfWestLongitudeDeg, dfSouthLatitudeDeg,
                         dfEastLongitudeDeg, dfNorthLatitudeDeg);
            }
            else
            {
                CPLDebug("GDAL", "Could not compute area of interest");
                dfWestLongitudeDeg = 0;
                dfSouthLatitudeDeg = 0;
                dfEastLongitudeDeg = 0;
                dfNorthLatitudeDeg = 0;
            }
            delete poCT;
        }

        delete poGeog;
    }

    return ret;
}

/************************************************************************/
/*                    GDALGCPAntimeridianUnwrap()                       */
/************************************************************************/

/* Deal with discontinuties of dfGCPX longitudes around the anti-meridian.
 * Cf https://github.com/OSGeo/gdal/issues/8371
 */
static void GDALGCPAntimeridianUnwrap(int nGCPCount, GDAL_GCP *pasGCPList,
                                      const OGRSpatialReference &oSRS,
                                      CSLConstList papszOptions)
{
    const char *pszGCPAntimeridianUnwrap =
        CSLFetchNameValueDef(papszOptions, "GCP_ANTIMERIDIAN_UNWRAP", "AUTO");
    const bool bForced = EQUAL(pszGCPAntimeridianUnwrap, "YES") ||
                         EQUAL(pszGCPAntimeridianUnwrap, "ON") ||
                         EQUAL(pszGCPAntimeridianUnwrap, "TRUE") ||
                         EQUAL(pszGCPAntimeridianUnwrap, "1");
    if (bForced || (!oSRS.IsEmpty() && oSRS.IsGeographic() &&
                    fabs(oSRS.GetAngularUnits(nullptr) -
                         CPLAtof(SRS_UA_DEGREE_CONV)) < 1e-8 &&
                    EQUAL(pszGCPAntimeridianUnwrap, "AUTO")))
    {
        if (!bForced)
        {
            // Proceed to unwrapping only if the longitudes are within
            // [-180, -170] or [170, 180]
            for (int i = 0; i < nGCPCount; ++i)
            {
                const double dfLongAbs = fabs(pasGCPList[i].dfGCPX);
                if (dfLongAbs > 180 || dfLongAbs < 170)
                {
                    return;
                }
            }
        }

        bool bDone = false;
        for (int i = 0; i < nGCPCount; ++i)
        {
            if (pasGCPList[i].dfGCPX < 0)
            {
                if (!bDone)
                {
                    bDone = true;
                    CPLDebug("WARP", "GCP longitude unwrapping");
                }
                pasGCPList[i].dfGCPX += 360;
            }
        }
    }
}

/************************************************************************/
/*              GDALGetGenImgProjTranformerOptionList()                 */
/************************************************************************/

/** Return a XML string describing options accepted by
 * GDALCreateGenImgProjTransformer2().
 *
 * @since 3.11
 */
const char *GDALGetGenImgProjTranformerOptionList(void)
{
    return "<OptionList>"
           "<Option name='SRC_SRS' type='string' description='WKT SRS, or any "
           "string recognized by OGRSpatialReference::SetFromUserInput(), to "
           "be used as an override for CRS of input dataset'/>"
           "<Option name='DST_SRS' type='string' description='WKT SRS, or any "
           "string recognized by OGRSpatialReference::SetFromUserInput(), to "
           "be used as an override for CRS of output dataset'/>"
           "<Option name='PROMOTE_TO_3D' type='boolean' description='"
           "Whether to promote SRC_SRS / DST_SRS to 3D.' "
           "default='NO'/>"
           "<Option name='COORDINATE_OPERATION' type='string' description='"
           "Coordinate operation, as a PROJ or WKT string, used as an override "
           "over the normally computed pipeline. The pipeline must take into "
           "account the axis order of the source and target SRS.'/>"
           "<Option name='ALLOW_BALLPARK' type='boolean' description='"
           "Whether ballpark coordinate operations are allowed.' "
           "default='YES'/>"
           "<Option name='ONLY_BEST' type='string-select' "
           "description='"
           "By default (at least in the PROJ 9.x series), PROJ may use "
           "coordinate operations that are not the \"best\" if resources "
           "(typically grids) needed to use them are missing. It will then "
           "fallback to other coordinate operations that have a lesser "
           "accuracy, for example using Helmert transformations, or in the "
           "absence of such operations, to ones with potential very rough "
           " accuracy, using \"ballpark\" transformations (see "
           "https://proj.org/glossary.html). "
           "When calling this method with YES, PROJ will only consider the "
           "\"best\" operation, and error out (at Transform() time) if they "
           "cannot be used. This method may be used together with "
           "ALLOW_BALLPARK=NO to only allow best operations that have a known "
           "accuracy. Note that this method has no effect on PROJ versions "
           "before 9.2. The default value for this option can be also set with "
           "the PROJ_ONLY_BEST_DEFAULT environment variable, or with the "
           "\"only_best_default\" setting of proj.ini. Setting "
           "ONLY_BEST=YES/NO overrides such default value' default='AUTO'>"
           "  <Value>AUTO</Value>"
           "  <Value>YES</Value>"
           "  <Value>NO</Value>"
           "</Option>"
           "<Option name='COORDINATE_EPOCH' type='float' description='"
           "Coordinate epoch, expressed as a decimal year. Useful for "
           "time-dependent coordinate operations.'/>"
           "<Option name='SRC_COORDINATE_EPOCH' type='float' description='"
           "Coordinate epoch of source CRS, expressed as a decimal year. "
           "Useful for time-dependent coordinate operations.'/>"
           "<Option name='DST_COORDINATE_EPOCH' type='float' description='"
           "Coordinate epoch of target CRS, expressed as a decimal year. "
           "Useful for time-dependent coordinate operations.'/>"
           "<Option name='GCPS_OK' type='boolean' description='"
           "Allow use of GCPs.' default='YES'/>"
           "<Option name='REFINE_MINIMUM_GCPS' type='int' description='"
           "The minimum amount of GCPs that should be available after the "
           "refinement'/>"
           "<Option name='REFINE_TOLERANCE' type='float' description='"
           "The tolerance that specifies when a GCP will be eliminated.'/>"
           "<Option name='MAX_GCP_ORDER' type='int' description='"
           "The maximum order to use for GCP derived polynomials if possible. "
           "The default is to autoselect based on the number of GCPs. A value "
           "of -1 triggers use of Thin Plate Spline instead of polynomials.'/>"
           "<Option name='GCP_ANTIMERIDIAN_UNWRAP' type='string-select' "
           "description='"
           "Whether to \"unwrap\" longitudes of ground control points that "
           "span the antimeridian. For datasets with GCPs in "
           "longitude/latitude coordinate space spanning the antimeridian, "
           "longitudes will have a discontinuity on +/- 180 deg, and will "
           "result in a subset of the GCPs with longitude in the [-180,-170] "
           "range and another subset in [170, 180]. By default (AUTO), that "
           "situation will be detected and longitudes in [-180,-170] will be "
           "shifted to [180, 190] to get a continuous set. This option can be "
           "set to YES to force that behavior (useful if no SRS information is "
           "available), or to NO to disable it.' default='AUTO'>"
           "  <Value>AUTO</Value>"
           "  <Value>YES</Value>"
           "  <Value>NO</Value>"
           "</Option>"
           "<Option name='SRC_METHOD' alias='METHOD' type='string-select' "
           "description='"
           "Force only one geolocation method to be considered on the source "
           "dataset. Will be used for pixel/line to georef transformation on "
           "the source dataset. NO_GEOTRANSFORM can be used to specify the "
           "identity geotransform (ungeoreferenced image)'>"
           "  <Value>GEOTRANSFORM</Value>"
           "  <Value>GCP_POLYNOMIAL</Value>"
           "  <Value>GCP_TPS</Value>"
           "  <Value>GCP_HOMOGRAPHY</Value>"
           "  <Value>GEOLOC_ARRAY</Value>"
           "  <Value>RPC</Value>"
           "  <Value>NO_GEOTRANSFORM</Value>"
           "</Option>"
           "<Option name='DST_METHOD' type='string-select' description='"
           "Force only one geolocation method to be considered on the target "
           "dataset. Will be used for pixel/line to georef transformation on "
           "the targe dataset. NO_GEOTRANSFORM can be used to specify the "
           "identity geotransform (ungeoreferenced image)'>"
           "  <Value>GEOTRANSFORM</Value>"
           "  <Value>GCP_POLYNOMIAL</Value>"
           "  <Value>GCP_TPS</Value>"
           "  <Value>GCP_HOMOGRAPHY</Value>"
           "  <Value>GEOLOC_ARRAY</Value>"
           "  <Value>RPC</Value>"
           "  <Value>NO_GEOTRANSFORM</Value>"
           "</Option>"
           "<Option name='RPC_HEIGHT' type='float' description='"
           "A fixed height to be used with RPC calculations. If RPC_HEIGHT and "
           "RPC_DEM are not specified but that the RPC metadata domain contains"
           " a HEIGHT_DEFAULT item (for example, the DIMAP driver may fill it),"
           "this value will be used as the RPC_HEIGHT. Otherwise, if none of "
           "RPC_HEIGHT and RPC_DEM are specified as transformer options and "
           "if HEIGHT_DEFAULT is no available, a height of 0 will be used.'/>"
           "<Option name='RPC_DEM' type='string' description='"
           "Name of a GDAL dataset (a DEM file typically) used to extract "
           "elevation offsets from. In this situation the Z passed into the "
           "transformation function is assumed to be height above ground. "
           "This option should be used in replacement of RPC_HEIGHT to provide "
           "a way of defining a non uniform ground for the target scene.'/>"
           "<Option name='RPC_HEIGHT_SCALE' type='float' description='"
           "Factor used to multiply heights above ground. Useful when "
           "elevation offsets of the DEM are not expressed in meters.'/>"
           "<Option name='RPC_DEMINTERPOLATION' type='string-select' "
           "description='DEM interpolation method' default='BILINEAR'>"
           "  <Value>NEAR</Value>"
           "  <Value>BILINEAR</Value>"
           "  <Value>CUBIC</Value>"
           "</Option>"
           "<Option name='RPC_DEM_MISSING_VALUE' type='float' description='"
           "Value of DEM height that must be used in case the DEM has nodata "
           "value at the sampling point, or if its extent does not cover the "
           "requested coordinate. When not specified, missing values will "
           "cause a failed transform.'/>"
           "<Option name='RPC_DEM_SRS' type='string' description='"
           "WKT SRS, or any string recognized by "
           "OGRSpatialReference::SetFromUserInput(), to be used as an "
           "override for DEM SRS. Useful if DEM SRS does not have an explicit "
           "vertical component.'/>"
           "<Option name='RPC_DEM_APPLY_VDATUM_SHIFT' type='boolean' "
           "description='"
           "Whether the vertical component of a compound SRS for the DEM "
           "should be used (when it is present). This is useful so as to "
           "be able to transform the raw values from the DEM expressed with "
           "respect to a geoid to the heights with respect to the WGS84 "
           "ellipsoid. When this is enabled, the GTIFF_REPORT_COMPD_CS "
           "configuration option will be also set temporarily so as to get "
           "the vertical information from GeoTIFF files.' default='YES'/>"
           "<Option name='RPC_PIXEL_ERROR_THRESHOLD' type='float' description='"
           "Overrides the dfPixErrThreshold parameter, i.e. the error "
           "(measured in pixels) allowed in the iterative solution of "
           "pixel/line to lat/long computations (the other way is always "
           "exact given the equations).'/>"
           "<Option name='RPC_MAX_ITERATIONS' type='int' description='"
           "Maximum number of iterations allowed in the iterative solution of "
           "pixel/line to lat/long computations. Default value is 10 in the "
           "absence of a DEM, or 20 if there is a DEM.'/>"
           "<Option name='RPC_FOOTPRINT' type='string' description='"
           "WKT or GeoJSON polygon (in long / lat coordinate space) with a "
           "validity footprint for the RPC. Any coordinate transformation that "
           "goes from or arrive outside this footprint will be considered "
           "invalid. This* is useful in situations where the RPC values become "
           "highly unstable outside of the area on which they have been "
           "computed for, potentially leading to undesirable \"echoes\" / "
           "false positives. This requires GDAL to be built against GEOS..'/>"
           "<Option name='RPC_MAX_ITERATIONS' type='int' description='"
           "Maximum number of iterations allowed in the iterative solution of "
           "pixel/line to lat/long computations. Default value is 10 in the "
           "absence of a DEM, or 20 if there is a DEM.'/>"
           "<Option name='INSERT_CENTER_LONG' type='boolean' description='"
           "May be set to FALSE to disable setting up a CENTER_LONG value on "
           "the coordinate system to rewrap things around the center of the "
           "image.' default='YES'/>"
           "<Option name='SRC_APPROX_ERROR_IN_SRS_UNIT' type='float' "
           "description='"
           "Use an approximate transformer for the source transformer. Must be "
           "defined together with SRC_APPROX_ERROR_IN_PIXEL to be taken into "
           "account.'/>"
           "<Option name='SRC_APPROX_ERROR_IN_PIXEL' type='float' "
           "description='"
           "Use an approximate transformer for the source transformer. Must be "
           "defined together with SRC_APPROX_ERROR_IN_SRS_UNIT to be taken "
           "into "
           "account.'/>"
           "<Option name='DST_APPROX_ERROR_IN_SRS_UNIT' type='float' "
           "description='"
           "Use an approximate transformer for the target transformer. Must be "
           "defined together with DST_APPROX_ERROR_IN_PIXEL to be taken into "
           "account.'/>"
           "<Option name='DST_APPROX_ERROR_IN_PIXEL' type='float' "
           "description='"
           "Use an approximate transformer for the target transformer. Must be "
           "defined together with DST_APPROX_ERROR_IN_SRS_UNIT to be taken "
           "into "
           "account.'/>"
           "<Option name='REPROJECTION_APPROX_ERROR_IN_SRC_SRS_UNIT' "
           "type='float' "
           "description='"
           "Use an approximate transformer for the coordinate reprojection. "
           "Must be used together with "
           "REPROJECTION_APPROX_ERROR_IN_DST_SRS_UNIT to be taken into "
           "account.'/>"
           "<Option name='REPROJECTION_APPROX_ERROR_IN_DST_SRS_UNIT' "
           "type='float' "
           "description='"
           "Use an approximate transformer for the coordinate reprojection. "
           "Must be used together with "
           "REPROJECTION_APPROX_ERROR_IN_SRC_SRS_UNIT to be taken into "
           "account.'/>"
           "<Option name='AREA_OF_INTEREST' type='string' "
           "description='"
           "Area of interest, as "
           "west_lon_deg,south_lat_deg,east_lon_deg,north_lat_deg, used to "
           "compute the best coordinate operation between the source and "
           "target SRS. If not specified, the bounding box of the source "
           "raster will be used.'/>"
           "<Option name='GEOLOC_BACKMAP_OVERSAMPLE_FACTOR' type='float' "
           "min='0.1' max='2' description='"
           "Oversample factor used to derive the size of the \"backmap\" used "
           "for geolocation array transformers.' default='1.3'/>"
           "<Option name='GEOLOC_USE_TEMP_DATASETS' type='boolean' "
           "description='"
           "Whether temporary GeoTIFF datasets should be used to store the "
           "backmap. The default is NO, that is to use in-memory arrays, "
           "unless the number of pixels of the geolocation array is greater "
           "than 16 megapixels.' default='NO'/>"
           "<Option name='GEOLOC_ARRAY' alias='SRC_GEOLOC_ARRAY' type='string' "
           "description='"
           "Name of a GDAL dataset containing a geolocation array and "
           "associated metadata. This is an alternative to having geolocation "
           "information described in the GEOLOCATION metadata domain of the "
           "source dataset. The dataset specified may have a GEOLOCATION "
           "metadata domain containing appropriate metadata, however default "
           "values are assigned for all omitted items. X_BAND defaults to 1 "
           "and Y_BAND to 2, however the dataset must contain exactly 2 bands. "
           "PIXEL_OFFSET and LINE_OFFSET default to 0. PIXEL_STEP and "
           "LINE_STEP default to the ratio of the width/height of the source "
           "dataset divided by the with/height of the geolocation array. "
           "SRS defaults to the spatial reference system of the geolocation "
           "array dataset, if set, otherwise WGS84 is used. "
           "GEOREFERENCING_CONVENTION is selected from the main metadata "
           "domain if it is omitted from the GEOLOCATION domain, and if not "
           "available TOP_LEFT_CORNER is assigned as a default. "
           "If GEOLOC_ARRAY is set SRC_METHOD defaults to GEOLOC_ARRAY.'/>"
           "<Option name='DST_GEOLOC_ARRAY' type='string' "
           "description='"
           "Name of a GDAL dataset that contains at least 2 bands with the X "
           "and Y geolocation bands. This is an alternative to having "
           "geolocation information described in the GEOLOCATION metadata "
           "domain of the destination dataset. See SRC_GEOLOC_ARRAY "
           "description for details, assumptions, and defaults. If this "
           "option is set, DST_METHOD=GEOLOC_ARRAY will be assumed if not "
           "set.'/>"
           "<Option name='GEOLOC_NORMALIZE_LONGITUDE_MINUS_180_PLUS_180' "
           "type='boolean' "
           "description='"
           "Force geolocation longitudes into -180,180 when longitude/latitude "
           "is the coordinate system of the geolocation arrays' default='NO'>"
           "  <Value>YES</Value>"
           "  <Value>NO</Value>"
           "</Option>"
           "<Option name='NUM_THREADS' type='string' "
           "description='Number of threads to use'/>"
           "</OptionList>";
}

/************************************************************************/
/*                  GDALCreateGenImgProjTransformer2()                  */
/************************************************************************/

/* clang-format off */
/**
 * Create image to image transformer.
 *
 * This function creates a transformation object that maps from pixel/line
 * coordinates on one image to pixel/line coordinates on another image.  The
 * images may potentially be georeferenced in different coordinate systems,
 * and may used GCPs to map between their pixel/line coordinates and
 * georeferenced coordinates (as opposed to the default assumption that their
 * geotransform should be used).
 *
 * This transformer potentially performs three concatenated transformations.
 *
 * The first stage is from source image pixel/line coordinates to source
 * image georeferenced coordinates, and may be done using the geotransform,
 * or if not defined using a polynomial model derived from GCPs.  If GCPs
 * are used this stage is accomplished using GDALGCPTransform().
 *
 * The second stage is to change projections from the source coordinate system
 * to the destination coordinate system, assuming they differ.  This is
 * accomplished internally using GDALReprojectionTransform().
 *
 * The third stage is converting from destination image georeferenced
 * coordinates to destination image coordinates.  This is done using the
 * destination image geotransform, or if not available, using a polynomial
 * model derived from GCPs. If GCPs are used this stage is accomplished using
 * GDALGCPTransform().  This stage is skipped if hDstDS is NULL when the
 * transformation is created.
 *
 * Supported Options (specified with the -to switch of gdalwarp for example):
 * <ul>
 * <li> SRC_SRS: WKT SRS, or any string recognized by
 * OGRSpatialReference::SetFromUserInput(), to be used as an override for
 * hSrcDS.</li>
 * <li> DST_SRS: WKT SRS, or any string recognized by
 * OGRSpatialReference::SetFromUserInput(),  to be used as an override for
 * hDstDS.
 * </li>
 * <li>PROMOTE_TO_3D=YES/NO: whether to promote SRC_SRS / DST_SRS to 3D.
 * Default is NO</li>
 * <li> COORDINATE_OPERATION: (GDAL &gt;= 3.0) Coordinate operation, as
 * a PROJ or WKT string, used as an override over the normally computed
 * pipeline. The pipeline must take into account the axis order of the source
 * and target SRS.
 * </li>
 * <li> ALLOW_BALLPARK=YES/NO: (GDAL &gt;= 3.11) Whether ballpark coordinate
 * operations are allowed. Defaults to YES.</li>
 * <li> ONLY_BEST=YES/NO/AUTO: (GDAL &gt;= 3.11) By default (at least in the
 * PROJ 9.x series), PROJ may use coordinate
 * operations that are not the "best" if resources (typically grids) needed
 * to use them are missing. It will then fallback to other coordinate operations
 * that have a lesser accuracy, for example using Helmert transformations,
 * or in the absence of such operations, to ones with potential very rought
 * accuracy, using "ballpark" transformations
 * (see https://proj.org/glossary.html).
 * When calling this method with YES, PROJ will only consider the
 * "best" operation, and error out (at Transform() time) if they cannot be
 * used.
 * This method may be used together with ALLOW_BALLPARK=NO to
 * only allow best operations that have a known accuracy.
 * Note that this method has no effect on PROJ versions before 9.2.
 * The default value for this option can be also set with the
 * PROJ_ONLY_BEST_DEFAULT environment variable, or with the "only_best_default"
 * setting of proj.ini. Calling SetOnlyBest() overrides such default value.</li>
 * <li> COORDINATE_EPOCH: (GDAL &gt;= 3.0) Coordinate epoch,
 * expressed as a decimal year. Useful for time-dependent coordinate operations.
 * </li>
 * <li> SRC_COORDINATE_EPOCH: (GDAL &gt;= 3.4) Coordinate epoch of source CRS,
 * expressed as a decimal year. Useful for time-dependent coordinate operations.
 * </li>
 * <li> DST_COORDINATE_EPOCH: (GDAL &gt;= 3.4) Coordinate epoch of target CRS,
 * expressed as a decimal year. Useful for time-dependent coordinate operations.
 * </li>
 * <li> GCPS_OK: If false, GCPs will not be used, default is TRUE.
 * </li>
 * <li> REFINE_MINIMUM_GCPS: The minimum amount of GCPs that should be available
 * after the refinement.
 * </li>
 * <li> REFINE_TOLERANCE: The tolerance that specifies when a GCP will be
 * eliminated.
 * </li>
 * <li> MAX_GCP_ORDER: the maximum order to use for GCP derived polynomials if
 * possible.  The default is to autoselect based on the number of GCPs.
 * A value of -1 triggers use of Thin Plate Spline instead of polynomials.
 * </li>
 * <li>GCP_ANTIMERIDIAN_UNWRAP=AUTO/YES/NO. (GDAL &gt;= 3.8) Whether to
 * "unwrap" longitudes of ground control points that span the antimeridian.
 * For datasets with GCPs in longitude/latitude coordinate space spanning the
 * antimeridian, longitudes will have a discontinuity on +/- 180 deg, and
 * will result in a subset of the GCPs with longitude in the [-180,-170] range
 * and another subset in [170, 180]. By default (AUTO), that situation will be
 * detected and longitudes in [-180,-170] will be shifted to [180, 190] to get
 * a continuous set. This option can be set to YES to force that behavior
 * (useful if no SRS information is available), or to NO to disable it.
 * </li>
 * <li> SRC_METHOD: may have a value which is one of GEOTRANSFORM, GCP_HOMOGRAPHY,
 * GCP_POLYNOMIAL, GCP_TPS, GEOLOC_ARRAY, RPC to force only one geolocation
 * method to be considered on the source dataset. Will be used for pixel/line
 * to georef transformation on the source dataset. NO_GEOTRANSFORM can be
 * used to specify the identity geotransform (ungeoreferenced image)
 * </li>
 * <li> DST_METHOD: may have a value which is one of GEOTRANSFORM,
 * GCP_POLYNOMIAL, GCP_HOMOGRAPHY, GCP_TPS, GEOLOC_ARRAY (added in 3.5), RPC to
 * force only one
 * geolocation method to be considered on the target dataset.  Will be used for
 * pixel/line to georef transformation on the destination dataset.
 * NO_GEOTRANSFORM can be used to specify the identity geotransform
 * (ungeoreferenced image)
 * </li>
 * <li> RPC_HEIGHT: A fixed height to be used with RPC
 * calculations. If RPC_HEIGHT and RPC_DEM are not specified but that the RPC
 * metadata domain contains a HEIGHT_DEFAULT item (for example, the DIMAP driver
 * may fill it), this value will be used as the RPC_HEIGHT. Otherwise, if none
 * of RPC_HEIGHT and RPC_DEM are specified as transformer
 * options and if HEIGHT_DEFAULT is no available, a height of 0 will be used.
 * </li>
 * <li> RPC_DEM: The name of a DEM file to be used with RPC
 * calculations. See GDALCreateRPCTransformerV2() for more details.
 * </li>
 * <li> Other RPC related options. See GDALCreateRPCTransformerV2()
 * </li>
 * <li>
 * INSERT_CENTER_LONG: May be set to FALSE to disable setting up a CENTER_LONG
 * value on the coordinate system to rewrap things around the center of the
 * image.
 * </li>
 * <li> SRC_APPROX_ERROR_IN_SRS_UNIT=err_threshold_in_SRS_units. (GDAL
 * &gt;= 2.2) Use an approximate transformer for the source transformer. Must be
 * defined together with SRC_APPROX_ERROR_IN_PIXEL to be taken into account.
 * </li>
 * <li> SRC_APPROX_ERROR_IN_PIXEL=err_threshold_in_pixel. (GDAL &gt;= 2.2) Use
 * an approximate transformer for the source transformer.. Must be defined
 * together with SRC_APPROX_ERROR_IN_SRS_UNIT to be taken into account.
 * </li>
 * <li>
 * DST_APPROX_ERROR_IN_SRS_UNIT=err_threshold_in_SRS_units. (GDAL &gt;= 2.2) Use
 * an approximate transformer for the destination transformer. Must be defined
 * together with DST_APPROX_ERROR_IN_PIXEL to be taken into account.
 * </li>
 * <li>
 * DST_APPROX_ERROR_IN_PIXEL=err_threshold_in_pixel. (GDAL &gt;= 2.2) Use an
 * approximate transformer for the destination transformer. Must be defined
 * together with DST_APPROX_ERROR_IN_SRS_UNIT to be taken into account.
 * </li>
 * <li>
 * REPROJECTION_APPROX_ERROR_IN_SRC_SRS_UNIT=err_threshold_in_src_SRS_units.
 * (GDAL &gt;= 2.2) Use an approximate transformer for the coordinate
 * reprojection. Must be used together with
 * REPROJECTION_APPROX_ERROR_IN_DST_SRS_UNIT to be taken into account.
 * </li>
 * <li>
 * REPROJECTION_APPROX_ERROR_IN_DST_SRS_UNIT=err_threshold_in_dst_SRS_units.
 * (GDAL &gt;= 2.2) Use an approximate transformer for the coordinate
 * reprojection. Must be used together with
 * REPROJECTION_APPROX_ERROR_IN_SRC_SRS_UNIT to be taken into account.
 * </li>
 * <li>
 * AREA_OF_INTEREST=west_lon_deg,south_lat_deg,east_lon_deg,north_lat_deg. (GDAL
 * &gt;= 3.0) Area of interest, used to compute the best coordinate operation
 * between the source and target SRS. If not specified, the bounding box of the
 * source raster will be used.
 * </li>
 * <li> GEOLOC_BACKMAP_OVERSAMPLE_FACTOR=[0.1,2]. (GDAL &gt;= 3.5) Oversample
 * factor used to derive the size of the "backmap" used for geolocation array
 * transformers. Default value is 1.3.
 * </li>
 * <li> GEOLOC_USE_TEMP_DATASETS=YES/NO.
 * (GDAL &gt;= 3.5) Whether temporary GeoTIFF datasets should be used to store
 * the backmap. The default is NO, that is to use in-memory arrays, unless the
 * number of pixels of the geolocation array is greater than 16 megapixels.
 * </li>
 * <li>
 * GEOLOC_ARRAY/SRC_GEOLOC_ARRAY=filename. (GDAL &gt;= 3.5.2) Name of a GDAL
 * dataset containing a geolocation array and associated metadata. This is an
 * alternative to having geolocation information described in the GEOLOCATION
 * metadata domain of the source dataset. The dataset specified may have a
 * GEOLOCATION metadata domain containing appropriate metadata, however default
 * values are assigned for all omitted items. X_BAND defaults to 1 and Y_BAND to
 * 2, however the dataset must contain exactly 2 bands. PIXEL_OFFSET and
 * LINE_OFFSET default to 0. PIXEL_STEP and LINE_STEP default to the ratio of
 * the width/height of the source dataset divided by the with/height of the
 * geolocation array. SRS defaults to the geolocation array dataset's spatial
 * reference system if set, otherwise WGS84 is used.
 * GEOREFERENCING_CONVENTION is selected from the main metadata domain if it
 * is omitted from the GEOLOCATION domain, and if not available
 * TOP_LEFT_CORNER is assigned as a default.
 * If GEOLOC_ARRAY is set SRC_METHOD
 * defaults to GEOLOC_ARRAY.
 * </li>
 * <li>DST_GEOLOC_ARRAY=filename. (GDAL &gt;= 3.5.2) Name of a
 * GDAL dataset that contains at least 2 bands with the X and Y geolocation
 * bands. This is an alternative to having geolocation information described in
 * the GEOLOCATION metadata domain of the destination dataset. See
 * SRC_GEOLOC_ARRAY description for details, assumptions, and defaults. If this
 * option is set, DST_METHOD=GEOLOC_ARRAY will be assumed if not set.
 * </li>
 * <li>GEOLOC_NORMALIZE_LONGITUDE_MINUS_180_PLUS_180=YES/NO. (GDAL &gt;= 3.12.0)
 * Whether to force geolocation longitudes into -180,180 when longitude/latitude is
 * the coordinate system of the geolocation arrays. The default is to enable this mode
 * when the values in the geolocation array are in the -180,180, otherwise NO.
 * </li>
 * </ul>
 *
 * The use case for the *_APPROX_ERROR_* options is when defining an approximate
 * transformer on top of the GenImgProjTransformer globally is not practical.
 * Such a use case is when the source dataset has RPC with a RPC DEM. In such
 * case we don't want to use the approximate transformer on the RPC
 * transformation, as the RPC DEM generally involves non-linearities that the
 * approximate transformer will not detect. In such case, we must a
 * non-approximated GenImgProjTransformer, but it might be worthwhile to use
 * approximate sub- transformers, for example on coordinate reprojection. For
 * example if warping from a source dataset with RPC to a destination dataset
 * with a UTM projection, since the inverse UTM transformation is rather costly.
 * In which case, one can use the REPROJECTION_APPROX_ERROR_IN_SRC_SRS_UNIT and
 * REPROJECTION_APPROX_ERROR_IN_DST_SRS_UNIT options.
 *
 * The list of supported options can also be programmatically obtained with
 * GDALGetGenImgProjTranformerOptionList().
 *
 * @param hSrcDS source dataset, or NULL.
 * @param hDstDS destination dataset (or NULL).
 * @param papszOptions NULL-terminated list of string options (or NULL).
 *
 * @return handle suitable for use GDALGenImgProjTransform(), and to be
 * deallocated with GDALDestroyGenImgProjTransformer() or NULL on failure.
 */
/* clang-format on */

void *GDALCreateGenImgProjTransformer2(GDALDatasetH hSrcDS, GDALDatasetH hDstDS,
                                       CSLConstList papszOptions)

{
    GDALValidateOptions(GDALGetGenImgProjTranformerOptionList(), papszOptions,
                        "option", "transformer options");

    double dfWestLongitudeDeg = 0.0;
    double dfSouthLatitudeDeg = 0.0;
    double dfEastLongitudeDeg = 0.0;
    double dfNorthLatitudeDeg = 0.0;
    bool bHasAreaOfInterest = false;
    if (const char *pszAreaOfInterest =
            CSLFetchNameValue(papszOptions, "AREA_OF_INTEREST"))
    {
        const CPLStringList aosTokens(
            CSLTokenizeString2(pszAreaOfInterest, ", ", 0));
        if (aosTokens.size() == 4)
        {
            dfWestLongitudeDeg = CPLAtof(aosTokens[0]);
            dfSouthLatitudeDeg = CPLAtof(aosTokens[1]);
            dfEastLongitudeDeg = CPLAtof(aosTokens[2]);
            dfNorthLatitudeDeg = CPLAtof(aosTokens[3]);
            bHasAreaOfInterest = true;
        }
    }

    const char *pszCO = CSLFetchNameValue(papszOptions, "COORDINATE_OPERATION");

    const auto SetAxisMapping =
        [papszOptions](OGRSpatialReference &oSRS, const char *pszPrefix)
    {
        const char *pszMapping = CSLFetchNameValue(
            papszOptions, std::string(pszPrefix)
                              .append("_DATA_AXIS_TO_SRS_AXIS_MAPPING")
                              .c_str());
        if (pszMapping)
        {
            CPLStringList aosTokens(CSLTokenizeString2(pszMapping, ",", 0));
            std::vector<int> anMapping;
            for (int i = 0; i < aosTokens.size(); ++i)
                anMapping.push_back(atoi(aosTokens[i]));
            oSRS.SetDataAxisToSRSAxisMapping(anMapping);
        }
        else
        {
            const char *pszStrategy = CSLFetchNameValueDef(
                papszOptions,
                std::string(pszPrefix).append("_AXIS_MAPPING_STRATEGY").c_str(),
                "TRADITIONAL_GIS_ORDER");
            if (EQUAL(pszStrategy, "TRADITIONAL_GIS_ORDER"))
                oSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
            else if (EQUAL(pszStrategy, "AUTHORITY_COMPLIANT"))
                oSRS.SetAxisMappingStrategy(OAMS_AUTHORITY_COMPLIANT);
            else
            {
                CPLError(CE_Warning, CPLE_AppDefined,
                         "Unrecognized value '%s' for %s", pszStrategy,
                         std::string(pszPrefix)
                             .append("_AXIS_MAPPING_STRATEGY")
                             .c_str());
                return false;
            }
        }
        return true;
    };

    /* -------------------------------------------------------------------- */
    /*      Initialize the transform info.                                  */
    /* -------------------------------------------------------------------- */
    GDALGenImgProjTransformInfo *psInfo =
        GDALCreateGenImgProjTransformerInternal();

    const auto DealWithForwardOrInverse =
        [bHasAreaOfInterest, &dfWestLongitudeDeg, &dfSouthLatitudeDeg,
         &dfEastLongitudeDeg, &dfNorthLatitudeDeg, pszCO, papszOptions,
         &SetAxisMapping](GDALGenImgProjTransformPart &part, GDALDatasetH hDS,
                          const char *pszPrefix, OGRSpatialReference &oSRS,
                          bool &bCanUseGeoTransform)
    {
        const int nOrder =
            atoi(CSLFetchNameValueDef(papszOptions, "MAX_GCP_ORDER", "0"));

        const bool bGCPUseOK =
            CPLTestBool(CSLFetchNameValueDef(papszOptions, "GCPS_OK", "YES"));
        const int nMinimumGcps = atoi(
            CSLFetchNameValueDef(papszOptions, "REFINE_MINIMUM_GCPS", "-1"));

        const char *pszRefineTolerance =
            CSLFetchNameValue(papszOptions, "REFINE_TOLERANCE");
        const bool bRefine = pszRefineTolerance != nullptr;
        const double dfTolerance =
            pszRefineTolerance ? CPLAtof(pszRefineTolerance) : 0.0;

        const std::string osSRSOptionName =
            std::string(pszPrefix).append("_SRS");
        const char *pszSRS =
            CSLFetchNameValue(papszOptions, osSRSOptionName.c_str());
        if (pszSRS)
        {
            if (pszSRS[0] != '\0' &&
                oSRS.SetFromUserInput(pszSRS) != OGRERR_NONE)
            {
                CPLError(CE_Failure, CPLE_AppDefined,
                         "Failed to import coordinate system `%s'.", pszSRS);
                return false;
            }
            if (!SetAxisMapping(oSRS, osSRSOptionName.c_str()))
                return false;
        }

        CSLConstList papszMD = nullptr;
        GDALRPCInfoV2 sRPCInfo;

        bCanUseGeoTransform = false;

        const char *pszMethod = CSLFetchNameValue(
            papszOptions, std::string(pszPrefix).append("_METHOD").c_str());
        if (!pszMethod && EQUAL(pszPrefix, "SRC"))
            pszMethod = CSLFetchNameValue(papszOptions, "METHOD");

        const char *pszGeolocArray = CSLFetchNameValue(
            papszOptions,
            std::string(pszPrefix).append("_GEOLOC_ARRAY").c_str());
        if (!pszGeolocArray && EQUAL(pszPrefix, "SRC"))
            pszGeolocArray = CSLFetchNameValue(papszOptions, "GEOLOC_ARRAY");
        if (!pszMethod && pszGeolocArray != nullptr)
            pszMethod = "GEOLOC_ARRAY";

        /* -------------------------------------------------------------------- */
        /*      Get forward and inverse geotransform for the source image.      */
        /* -------------------------------------------------------------------- */
        if (hDS == nullptr ||
            (pszMethod != nullptr && EQUAL(pszMethod, "NO_GEOTRANSFORM")))
        {
            part.adfGeoTransform[0] = 0.0;
            part.adfGeoTransform[1] = 1.0;
            part.adfGeoTransform[2] = 0.0;
            part.adfGeoTransform[3] = 0.0;
            part.adfGeoTransform[4] = 0.0;
            part.adfGeoTransform[5] = 1.0;
            memcpy(part.adfInvGeoTransform, part.adfGeoTransform,
                   sizeof(double) * 6);
        }
        else if ((pszMethod == nullptr || EQUAL(pszMethod, "GEOTRANSFORM")) &&
                 GDALGetGeoTransform(hDS, part.adfGeoTransform) == CE_None)
        {
            if (!GDALInvGeoTransform(part.adfGeoTransform,
                                     part.adfInvGeoTransform))
            {
                CPLError(CE_Failure, CPLE_AppDefined,
                         "Cannot invert geotransform");
                return false;
            }
            if (pszSRS == nullptr)
            {
                auto hSRS = GDALGetSpatialRef(hDS);
                if (hSRS)
                    oSRS = *(OGRSpatialReference::FromHandle(hSRS));
            }
            if (EQUAL(pszPrefix, "SRC"))
            {
                if (!bHasAreaOfInterest && pszCO == nullptr && !oSRS.IsEmpty())
                {
                    GDALComputeAreaOfInterest(
                        &oSRS, part.adfGeoTransform, GDALGetRasterXSize(hDS),
                        GDALGetRasterYSize(hDS), dfWestLongitudeDeg,
                        dfSouthLatitudeDeg, dfEastLongitudeDeg,
                        dfNorthLatitudeDeg);
                }
                bCanUseGeoTransform = true;
            }
        }
        else if (bGCPUseOK &&
                 ((pszMethod == nullptr && GDALGetGCPCount(hDS) >= 4 &&
                   GDALGetGCPCount(hDS) < 6) ||
                  (pszMethod != nullptr &&
                   EQUAL(pszMethod, "GCP_HOMOGRAPHY"))) &&
                 GDALGetGCPCount(hDS) > 0)
        {
            if (pszSRS == nullptr)
            {
                auto hSRS = GDALGetGCPSpatialRef(hDS);
                if (hSRS)
                    oSRS = *(OGRSpatialReference::FromHandle(hSRS));
            }

            const auto nGCPCount = GDALGetGCPCount(hDS);
            auto pasGCPList = GDALDuplicateGCPs(nGCPCount, GDALGetGCPs(hDS));
            GDALGCPAntimeridianUnwrap(nGCPCount, pasGCPList, oSRS,
                                      papszOptions);

            part.pTransformArg =
                GDALCreateHomographyTransformerFromGCPs(nGCPCount, pasGCPList);

            GDALDeinitGCPs(nGCPCount, pasGCPList);
            CPLFree(pasGCPList);

            if (part.pTransformArg == nullptr)
            {
                return false;
            }
            part.pTransformer = GDALHomographyTransform;
        }
        else if (bGCPUseOK &&
                 (pszMethod == nullptr || EQUAL(pszMethod, "GCP_POLYNOMIAL")) &&
                 GDALGetGCPCount(hDS) > 0 && nOrder >= 0)
        {
            if (pszSRS == nullptr)
            {
                auto hSRS = GDALGetGCPSpatialRef(hDS);
                if (hSRS)
                    oSRS = *(OGRSpatialReference::FromHandle(hSRS));
            }

            const auto nGCPCount = GDALGetGCPCount(hDS);
            auto pasGCPList = GDALDuplicateGCPs(nGCPCount, GDALGetGCPs(hDS));
            GDALGCPAntimeridianUnwrap(nGCPCount, pasGCPList, oSRS,
                                      papszOptions);

            if (bRefine)
            {
                part.pTransformArg = GDALCreateGCPRefineTransformer(
                    nGCPCount, pasGCPList, nOrder, FALSE, dfTolerance,
                    nMinimumGcps);
            }
            else
            {
                part.pTransformArg = GDALCreateGCPTransformer(
                    nGCPCount, pasGCPList, nOrder, FALSE);
            }

            GDALDeinitGCPs(nGCPCount, pasGCPList);
            CPLFree(pasGCPList);

            if (part.pTransformArg == nullptr)
            {
                return false;
            }
            part.pTransformer = GDALGCPTransform;
        }

        else if (bGCPUseOK && GDALGetGCPCount(hDS) > 0 && nOrder <= 0 &&
                 (pszMethod == nullptr || EQUAL(pszMethod, "GCP_TPS")))
        {
            if (pszSRS == nullptr)
            {
                auto hSRS = GDALGetGCPSpatialRef(hDS);
                if (hSRS)
                    oSRS = *(OGRSpatialReference::FromHandle(hSRS));
            }

            const auto nGCPCount = GDALGetGCPCount(hDS);
            auto pasGCPList = GDALDuplicateGCPs(nGCPCount, GDALGetGCPs(hDS));
            GDALGCPAntimeridianUnwrap(nGCPCount, pasGCPList, oSRS,
                                      papszOptions);

            part.pTransformArg = GDALCreateTPSTransformerInt(
                nGCPCount, pasGCPList, FALSE, papszOptions);

            GDALDeinitGCPs(nGCPCount, pasGCPList);
            CPLFree(pasGCPList);

            if (part.pTransformArg == nullptr)
            {
                return false;
            }
            part.pTransformer = GDALTPSTransform;
        }

        else if ((pszMethod == nullptr || EQUAL(pszMethod, "RPC")) &&
                 (papszMD = GDALGetMetadata(hDS, "RPC")) != nullptr &&
                 GDALExtractRPCInfoV2(papszMD, &sRPCInfo))
        {
            CPLStringList aosOptions(papszOptions);
            if (!CSLFetchNameValue(papszOptions, "RPC_HEIGHT") &&
                !CSLFetchNameValue(papszOptions, "RPC_DEM"))
            {
                if (const char *pszHEIGHT_DEFAULT =
                        CSLFetchNameValue(papszMD, "HEIGHT_DEFAULT"))
                {
                    CPLDebug("GDAL",
                             "For %s, using RPC_HEIGHT = HEIGHT_DEFAULT = %s",
                             pszPrefix, pszHEIGHT_DEFAULT);
                    aosOptions.SetNameValue("RPC_HEIGHT", pszHEIGHT_DEFAULT);
                }
            }
            part.pTransformArg = GDALCreateRPCTransformerV2(&sRPCInfo, FALSE, 0,
                                                            aosOptions.List());
            if (part.pTransformArg == nullptr)
            {
                return false;
            }
            part.pTransformer = GDALRPCTransform;
            if (pszSRS == nullptr)
            {
                oSRS.SetFromUserInput(SRS_WKT_WGS84_LAT_LONG);
                oSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
            }
        }

        else if ((pszMethod == nullptr || EQUAL(pszMethod, "GEOLOC_ARRAY")) &&
                 ((papszMD = GDALGetMetadata(hDS, "GEOLOCATION")) != nullptr ||
                  pszGeolocArray != nullptr))
        {
            CPLStringList aosGeolocMD;  // keep in this scope
            if (pszGeolocArray != nullptr)
            {
                if (papszMD != nullptr)
                {
                    CPLError(
                        CE_Warning, CPLE_AppDefined,
                        "Both GEOLOCATION metadata domain on the source "
                        "dataset "
                        "and [%s_]GEOLOC_ARRAY transformer option are set. "
                        "Only using the later.",
                        pszPrefix);
                }
                aosGeolocMD = GDALCreateGeolocationMetadata(
                    hDS, pszGeolocArray,
                    /* bIsSource= */ EQUAL(pszPrefix, "SRC"));
                if (aosGeolocMD.empty())
                {
                    return false;
                }
                papszMD = aosGeolocMD.List();
            }

            part.pTransformArg = GDALCreateGeoLocTransformerEx(
                hDS, papszMD, FALSE, nullptr, papszOptions);
            if (part.pTransformArg == nullptr)
            {
                return false;
            }
            part.pTransformer = GDALGeoLocTransform;
            if (pszSRS == nullptr)
            {
                pszSRS = CSLFetchNameValue(papszMD, "SRS");
                if (pszSRS)
                {
                    oSRS.SetFromUserInput(pszSRS);
                    oSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
                }
            }
        }

        else if (pszMethod != nullptr && EQUAL(pszPrefix, "SRC"))
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Unable to compute a %s based transformation between "
                     "pixel/line and georeferenced coordinates for %s.",
                     pszMethod, GDALGetDescription(hDS));

            return false;
        }

        else
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Unable to compute a transformation between pixel/line "
                     "and georeferenced coordinates for %s. "
                     "There is no affine transformation and no GCPs. "
                     "Specify transformation option %s_METHOD=NO_GEOTRANSFORM "
                     "to bypass this check.",
                     GDALGetDescription(hDS), pszPrefix);

            return false;
        }

        /* ---------------------------------------------------------------- */
        /*      Handle optional source approximation transformer.           */
        /* ---------------------------------------------------------------- */
        if (part.pTransformer)
        {
            const char *pszApproxErrorFwd = CSLFetchNameValue(
                papszOptions, std::string(pszPrefix)
                                  .append("_APPROX_ERROR_IN_SRS_UNIT")
                                  .c_str());
            const char *pszApproxErrorReverse = CSLFetchNameValue(
                papszOptions, std::string(pszPrefix)
                                  .append("_APPROX_ERROR_IN_PIXEL")
                                  .c_str());
            if (pszApproxErrorFwd && pszApproxErrorReverse)
            {
                void *pArg = GDALCreateApproxTransformer2(
                    part.pTransformer, part.pTransformArg,
                    CPLAtof(pszApproxErrorFwd), CPLAtof(pszApproxErrorReverse));
                if (pArg == nullptr)
                {
                    return false;
                }
                part.pTransformArg = pArg;
                part.pTransformer = GDALApproxTransform;
                GDALApproxTransformerOwnsSubtransformer(part.pTransformArg,
                                                        TRUE);
            }
        }

        return true;
    };

    /* -------------------------------------------------------------------- */
    /*      Get forward and inverse geotransform for the source image.      */
    /* -------------------------------------------------------------------- */
    bool bCanUseSrcGeoTransform = false;
    OGRSpatialReference oSrcSRS;
    if (!DealWithForwardOrInverse(psInfo->sSrcParams, hSrcDS, "SRC", oSrcSRS,
                                  bCanUseSrcGeoTransform))
    {
        GDALDestroyGenImgProjTransformer(psInfo);
        return nullptr;
    }

    /* -------------------------------------------------------------------- */
    /*      Get forward and inverse geotransform for destination image.     */
    /*      If we have no destination use a unit transform.                 */
    /* -------------------------------------------------------------------- */
    bool bIgnored = false;
    OGRSpatialReference oDstSRS;
    if (!DealWithForwardOrInverse(psInfo->sDstParams, hDstDS, "DST", oDstSRS,
                                  bIgnored))
    {
        GDALDestroyGenImgProjTransformer(psInfo);
        return nullptr;
    }

    /* -------------------------------------------------------------------- */
    /*      Setup reprojection.                                             */
    /* -------------------------------------------------------------------- */

    if (CPLFetchBool(papszOptions, "STRIP_VERT_CS", false))
    {
        if (oSrcSRS.IsCompound())
        {
            oSrcSRS.StripVertical();
        }
        if (oDstSRS.IsCompound())
        {
            oDstSRS.StripVertical();
        }
    }

    const bool bMayInsertCenterLong =
        (bCanUseSrcGeoTransform && !oSrcSRS.IsEmpty() && hSrcDS &&
         CPLFetchBool(papszOptions, "INSERT_CENTER_LONG", true));
    const char *pszSrcCoordEpoch =
        CSLFetchNameValue(papszOptions, "SRC_COORDINATE_EPOCH");
    const char *pszDstCoordEpoch =
        CSLFetchNameValue(papszOptions, "DST_COORDINATE_EPOCH");
    if ((!oSrcSRS.IsEmpty() && !oDstSRS.IsEmpty() &&
         (pszSrcCoordEpoch || pszDstCoordEpoch || !oSrcSRS.IsSame(&oDstSRS) ||
          (oSrcSRS.IsGeographic() && bMayInsertCenterLong))) ||
        pszCO)
    {
        CPLStringList aosOptions;

        if (bMayInsertCenterLong)
        {
            InsertCenterLong(hSrcDS, &oSrcSRS, aosOptions);
        }

        if (CPLFetchBool(papszOptions, "PROMOTE_TO_3D", false))
        {
            oSrcSRS.PromoteTo3D(nullptr);
            oDstSRS.PromoteTo3D(nullptr);
        }

        if (!(dfWestLongitudeDeg == 0.0 && dfSouthLatitudeDeg == 0.0 &&
              dfEastLongitudeDeg == 0.0 && dfNorthLatitudeDeg == 0.0))
        {
            aosOptions.SetNameValue(
                "AREA_OF_INTEREST",
                CPLSPrintf("%.16g,%.16g,%.16g,%.16g", dfWestLongitudeDeg,
                           dfSouthLatitudeDeg, dfEastLongitudeDeg,
                           dfNorthLatitudeDeg));
        }
        if (pszCO)
        {
            aosOptions.SetNameValue("COORDINATE_OPERATION", pszCO);
        }

        const char *pszCoordEpoch =
            CSLFetchNameValue(papszOptions, "COORDINATE_EPOCH");
        if (pszCoordEpoch)
        {
            aosOptions.SetNameValue("COORDINATE_EPOCH", pszCoordEpoch);
        }

        if (pszSrcCoordEpoch)
        {
            aosOptions.SetNameValue("SRC_COORDINATE_EPOCH", pszSrcCoordEpoch);
            oSrcSRS.SetCoordinateEpoch(CPLAtof(pszSrcCoordEpoch));
        }

        if (pszDstCoordEpoch)
        {
            aosOptions.SetNameValue("DST_COORDINATE_EPOCH", pszDstCoordEpoch);
            oDstSRS.SetCoordinateEpoch(CPLAtof(pszDstCoordEpoch));
        }

        if (const char *pszAllowBallpark =
                CSLFetchNameValue(papszOptions, "ALLOW_BALLPARK"))
        {
            aosOptions.SetNameValue("ALLOW_BALLPARK", pszAllowBallpark);
        }

        if (const char *pszOnlyBest =
                CSLFetchNameValue(papszOptions, "ONLY_BEST"))
        {
            aosOptions.SetNameValue("ONLY_BEST", pszOnlyBest);
        }

        psInfo->pReprojectArg = GDALCreateReprojectionTransformerEx(
            !oSrcSRS.IsEmpty() ? OGRSpatialReference::ToHandle(&oSrcSRS)
                               : nullptr,
            !oDstSRS.IsEmpty() ? OGRSpatialReference::ToHandle(&oDstSRS)
                               : nullptr,
            aosOptions.List());

        if (pszCO)
        {
            psInfo->bHasCustomTransformationPipeline = true;
        }

        if (psInfo->pReprojectArg == nullptr)
        {
            GDALDestroyGenImgProjTransformer(psInfo);
            return nullptr;
        }
        psInfo->pReproject = GDALReprojectionTransform;

        /* --------------------------------------------------------------------
         */
        /*      Handle optional reprojection approximation transformer. */
        /* --------------------------------------------------------------------
         */
        const char *psApproxErrorFwd = CSLFetchNameValue(
            papszOptions, "REPROJECTION_APPROX_ERROR_IN_DST_SRS_UNIT");
        const char *psApproxErrorReverse = CSLFetchNameValue(
            papszOptions, "REPROJECTION_APPROX_ERROR_IN_SRC_SRS_UNIT");
        if (psApproxErrorFwd && psApproxErrorReverse)
        {
            void *pArg = GDALCreateApproxTransformer2(
                psInfo->pReproject, psInfo->pReprojectArg,
                CPLAtof(psApproxErrorFwd), CPLAtof(psApproxErrorReverse));
            if (pArg == nullptr)
            {
                GDALDestroyGenImgProjTransformer(psInfo);
                return nullptr;
            }
            psInfo->pReprojectArg = pArg;
            psInfo->pReproject = GDALApproxTransform;
            GDALApproxTransformerOwnsSubtransformer(psInfo->pReprojectArg,
                                                    TRUE);
        }
    }

    return psInfo;
}

/************************************************************************/
/*                  GDALRefreshGenImgProjTransformer()                  */
/************************************************************************/

void GDALRefreshGenImgProjTransformer(void *hTransformArg)
{
    GDALGenImgProjTransformInfo *psInfo =
        static_cast<GDALGenImgProjTransformInfo *>(hTransformArg);

    if (psInfo->pReprojectArg &&
        psInfo->bCheckWithInvertPROJ != GetCurrentCheckWithInvertPROJ())
    {
        psInfo->bCheckWithInvertPROJ = !psInfo->bCheckWithInvertPROJ;

        CPLXMLNode *psXML =
            GDALSerializeTransformer(psInfo->pReproject, psInfo->pReprojectArg);
        GDALDestroyTransformer(psInfo->pReprojectArg);
        GDALDeserializeTransformer(psXML, &psInfo->pReproject,
                                   &psInfo->pReprojectArg);
        CPLDestroyXMLNode(psXML);
    }
}

/************************************************************************/
/*                  GDALCreateGenImgProjTransformer3()                  */
/************************************************************************/

/**
 * Create image to image transformer.
 *
 * This function creates a transformation object that maps from pixel/line
 * coordinates on one image to pixel/line coordinates on another image.  The
 * images may potentially be georeferenced in different coordinate systems,
 * and may used GCPs to map between their pixel/line coordinates and
 * georeferenced coordinates (as opposed to the default assumption that their
 * geotransform should be used).
 *
 * This transformer potentially performs three concatenated transformations.
 *
 * The first stage is from source image pixel/line coordinates to source
 * image georeferenced coordinates, and may be done using the geotransform,
 * or if not defined using a polynomial model derived from GCPs.  If GCPs
 * are used this stage is accomplished using GDALGCPTransform().
 *
 * The second stage is to change projections from the source coordinate system
 * to the destination coordinate system, assuming they differ.  This is
 * accomplished internally using GDALReprojectionTransform().
 *
 * The third stage is converting from destination image georeferenced
 * coordinates to destination image coordinates.  This is done using the
 * destination image geotransform, or if not available, using a polynomial
 * model derived from GCPs. If GCPs are used this stage is accomplished using
 * GDALGCPTransform().  This stage is skipped if hDstDS is NULL when the
 * transformation is created.
 *
 * @param pszSrcWKT source WKT (or NULL).
 * @param padfSrcGeoTransform source geotransform (or NULL).
 * @param pszDstWKT destination WKT (or NULL).
 * @param padfDstGeoTransform destination geotransform (or NULL).
 *
 * @return handle suitable for use GDALGenImgProjTransform(), and to be
 * deallocated with GDALDestroyGenImgProjTransformer() or NULL on failure.
 */

void *GDALCreateGenImgProjTransformer3(const char *pszSrcWKT,
                                       const double *padfSrcGeoTransform,
                                       const char *pszDstWKT,
                                       const double *padfDstGeoTransform)

{
    OGRSpatialReference oSrcSRS;
    if (pszSrcWKT)
    {
        oSrcSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        if (pszSrcWKT[0] != '\0' &&
            oSrcSRS.importFromWkt(pszSrcWKT) != OGRERR_NONE)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Failed to import coordinate system `%s'.", pszSrcWKT);
            return nullptr;
        }
    }

    OGRSpatialReference oDstSRS;
    if (pszDstWKT)
    {
        oDstSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        if (pszDstWKT[0] != '\0' &&
            oDstSRS.importFromWkt(pszDstWKT) != OGRERR_NONE)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Failed to import coordinate system `%s'.", pszDstWKT);
            return nullptr;
        }
    }
    return GDALCreateGenImgProjTransformer4(
        OGRSpatialReference::ToHandle(&oSrcSRS), padfSrcGeoTransform,
        OGRSpatialReference::ToHandle(&oDstSRS), padfDstGeoTransform, nullptr);
}

/************************************************************************/
/*                  GDALCreateGenImgProjTransformer4()                  */
/************************************************************************/

/**
 * Create image to image transformer.
 *
 * Similar to GDALCreateGenImgProjTransformer3(), except that it takes
 * OGRSpatialReferenceH objects and options.
 * The options are the ones supported by GDALCreateReprojectionTransformerEx()
 *
 * @since GDAL 3.0
 */
void *GDALCreateGenImgProjTransformer4(OGRSpatialReferenceH hSrcSRS,
                                       const double *padfSrcGeoTransform,
                                       OGRSpatialReferenceH hDstSRS,
                                       const double *padfDstGeoTransform,
                                       const char *const *papszOptions)
{
    /* -------------------------------------------------------------------- */
    /*      Initialize the transform info.                                  */
    /* -------------------------------------------------------------------- */
    GDALGenImgProjTransformInfo *psInfo =
        GDALCreateGenImgProjTransformerInternal();

    /* -------------------------------------------------------------------- */
    /*      Get forward and inverse geotransform for the source image.      */
    /* -------------------------------------------------------------------- */

    const auto SetParams =
        [](GDALGenImgProjTransformPart &part, const double *padfGT)
    {
        if (padfGT)
        {
            memcpy(part.adfGeoTransform, padfGT, sizeof(part.adfGeoTransform));
            if (!GDALInvGeoTransform(part.adfGeoTransform,
                                     part.adfInvGeoTransform))
            {
                CPLError(CE_Failure, CPLE_AppDefined,
                         "Cannot invert geotransform");
                return false;
            }
        }
        else
        {
            part.adfGeoTransform[0] = 0.0;
            part.adfGeoTransform[1] = 1.0;
            part.adfGeoTransform[2] = 0.0;
            part.adfGeoTransform[3] = 0.0;
            part.adfGeoTransform[4] = 0.0;
            part.adfGeoTransform[5] = 1.0;
            memcpy(part.adfInvGeoTransform, part.adfGeoTransform,
                   sizeof(double) * 6);
        }
        return true;
    };

    if (!SetParams(psInfo->sSrcParams, padfSrcGeoTransform))
    {
        GDALDestroyGenImgProjTransformer(psInfo);
        return nullptr;
    }

    /* -------------------------------------------------------------------- */
    /*      Setup reprojection.                                             */
    /* -------------------------------------------------------------------- */
    OGRSpatialReference *poSrcSRS = OGRSpatialReference::FromHandle(hSrcSRS);
    OGRSpatialReference *poDstSRS = OGRSpatialReference::FromHandle(hDstSRS);
    if (!poSrcSRS->IsEmpty() && !poDstSRS->IsEmpty() &&
        !poSrcSRS->IsSame(poDstSRS))
    {
        psInfo->pReprojectArg =
            GDALCreateReprojectionTransformerEx(hSrcSRS, hDstSRS, papszOptions);
        if (psInfo->pReprojectArg == nullptr)
        {
            GDALDestroyGenImgProjTransformer(psInfo);
            return nullptr;
        }
        psInfo->pReproject = GDALReprojectionTransform;
    }

    /* -------------------------------------------------------------------- */
    /*      Get forward and inverse geotransform for destination image.     */
    /*      If we have no destination matrix use a unit transform.          */
    /* -------------------------------------------------------------------- */
    if (!SetParams(psInfo->sDstParams, padfDstGeoTransform))
    {
        GDALDestroyGenImgProjTransformer(psInfo);
        return nullptr;
    }

    return psInfo;
}

/************************************************************************/
/*            GDALSetGenImgProjTransformerDstGeoTransform()             */
/************************************************************************/

/**
 * Set GenImgProj output geotransform.
 *
 * Normally the "destination geotransform", or transformation between
 * georeferenced output coordinates and pixel/line coordinates on the
 * destination file is extracted from the destination file by
 * GDALCreateGenImgProjTransformer() and stored in the GenImgProj private
 * info.  However, sometimes it is inconvenient to have an output file
 * handle with appropriate geotransform information when creating the
 * transformation.  For these cases, this function can be used to apply
 * the destination geotransform.
 *
 * @param hTransformArg the handle to update.
 * @param padfGeoTransform the destination geotransform to apply (six doubles).
 */

void GDALSetGenImgProjTransformerDstGeoTransform(void *hTransformArg,
                                                 const double *padfGeoTransform)

{
    VALIDATE_POINTER0(hTransformArg,
                      "GDALSetGenImgProjTransformerDstGeoTransform");

    GDALGenImgProjTransformInfo *psInfo =
        static_cast<GDALGenImgProjTransformInfo *>(hTransformArg);

    memcpy(psInfo->sDstParams.adfGeoTransform, padfGeoTransform,
           sizeof(double) * 6);
    if (!GDALInvGeoTransform(psInfo->sDstParams.adfGeoTransform,
                             psInfo->sDstParams.adfInvGeoTransform))
    {
        CPLError(CE_Failure, CPLE_AppDefined, "Cannot invert geotransform");
    }
}

/************************************************************************/
/*                  GDALDestroyGenImgProjTransformer()                  */
/************************************************************************/

/**
 * GenImgProjTransformer deallocator.
 *
 * This function is used to deallocate the handle created with
 * GDALCreateGenImgProjTransformer().
 *
 * @param hTransformArg the handle to deallocate.
 */

void GDALDestroyGenImgProjTransformer(void *hTransformArg)

{
    if (hTransformArg == nullptr)
        return;

    GDALGenImgProjTransformInfo *psInfo =
        static_cast<GDALGenImgProjTransformInfo *>(hTransformArg);

    if (psInfo->sSrcParams.pTransformArg != nullptr)
        GDALDestroyTransformer(psInfo->sSrcParams.pTransformArg);

    if (psInfo->sDstParams.pTransformArg != nullptr)
        GDALDestroyTransformer(psInfo->sDstParams.pTransformArg);

    if (psInfo->pReprojectArg != nullptr)
        GDALDestroyTransformer(psInfo->pReprojectArg);

    CPLFree(psInfo);
}

/************************************************************************/
/*                      GDALGenImgProjTransform()                       */
/************************************************************************/

/**
 * Perform general image reprojection transformation.
 *
 * Actually performs the transformation setup in
 * GDALCreateGenImgProjTransformer().  This function matches the signature
 * required by the GDALTransformerFunc(), and more details on the arguments
 * can be found in that topic.
 */

#ifdef DEBUG_APPROX_TRANSFORMER
int countGDALGenImgProjTransform = 0;
#endif

int GDALGenImgProjTransform(void *pTransformArgIn, int bDstToSrc,
                            int nPointCount, double *padfX, double *padfY,
                            double *padfZ, int *panSuccess)
{
    GDALGenImgProjTransformInfo *psInfo =
        static_cast<GDALGenImgProjTransformInfo *>(pTransformArgIn);

#ifdef DEBUG_APPROX_TRANSFORMER
    CPLAssert(nPointCount > 0);
    countGDALGenImgProjTransform += nPointCount;
#endif

    for (int i = 0; i < nPointCount; i++)
    {
        panSuccess[i] = (padfX[i] != HUGE_VAL && padfY[i] != HUGE_VAL);
    }

    int ret = TRUE;

    /* -------------------------------------------------------------------- */
    /*      Convert from src (dst) pixel/line to src (dst)                  */
    /*      georeferenced coordinates.                                      */
    /* -------------------------------------------------------------------- */
    {
        const auto params = bDstToSrc ? psInfo->sDstParams : psInfo->sSrcParams;
        const double *padfGeoTransform = params.adfGeoTransform;
        void *pTransformArg = params.pTransformArg;
        GDALTransformerFunc pTransformer = params.pTransformer;

        if (pTransformArg != nullptr)
        {
            if (!pTransformer(pTransformArg, FALSE, nPointCount, padfX, padfY,
                              padfZ, panSuccess))
                ret = FALSE;
        }
        else
        {
            for (int i = 0; i < nPointCount; i++)
            {
                if (!panSuccess[i])
                    continue;

                const double dfNewX = padfGeoTransform[0] +
                                      padfX[i] * padfGeoTransform[1] +
                                      padfY[i] * padfGeoTransform[2];
                const double dfNewY = padfGeoTransform[3] +
                                      padfX[i] * padfGeoTransform[4] +
                                      padfY[i] * padfGeoTransform[5];

                padfX[i] = dfNewX;
                padfY[i] = dfNewY;
            }
        }
    }

    /* -------------------------------------------------------------------- */
    /*      Reproject if needed.                                            */
    /* -------------------------------------------------------------------- */
    if (psInfo->pReprojectArg)
    {
        if (!psInfo->pReproject(psInfo->pReprojectArg, bDstToSrc, nPointCount,
                                padfX, padfY, padfZ, panSuccess))
            ret = FALSE;
    }

    /* -------------------------------------------------------------------- */
    /*      Convert dst (src) georef coordinates back to pixel/line.        */
    /* -------------------------------------------------------------------- */
    {
        const auto params = bDstToSrc ? psInfo->sSrcParams : psInfo->sDstParams;
        const double *padfInvGeoTransform = params.adfInvGeoTransform;
        void *pTransformArg = params.pTransformArg;
        GDALTransformerFunc pTransformer = params.pTransformer;

        if (pTransformArg != nullptr)
        {
            if (!pTransformer(pTransformArg, TRUE, nPointCount, padfX, padfY,
                              padfZ, panSuccess))
                ret = FALSE;
        }
        else
        {
            for (int i = 0; i < nPointCount; i++)
            {
                if (!panSuccess[i])
                    continue;

                const double dfNewX = padfInvGeoTransform[0] +
                                      padfX[i] * padfInvGeoTransform[1] +
                                      padfY[i] * padfInvGeoTransform[2];
                const double dfNewY = padfInvGeoTransform[3] +
                                      padfX[i] * padfInvGeoTransform[4] +
                                      padfY[i] * padfInvGeoTransform[5];

                padfX[i] = dfNewX;
                padfY[i] = dfNewY;
            }
        }
    }

    return ret;
}

/************************************************************************/
/*              GDALTransformLonLatToDestGenImgProjTransformer()        */
/************************************************************************/

int GDALTransformLonLatToDestGenImgProjTransformer(void *hTransformArg,
                                                   double *pdfX, double *pdfY)
{
    GDALGenImgProjTransformInfo *psInfo =
        static_cast<GDALGenImgProjTransformInfo *>(hTransformArg);

    if (psInfo->pReprojectArg == nullptr ||
        psInfo->pReproject != GDALReprojectionTransform)
        return false;

    GDALReprojectionTransformInfo *psReprojInfo =
        static_cast<GDALReprojectionTransformInfo *>(psInfo->pReprojectArg);
    if (psReprojInfo->poForwardTransform == nullptr ||
        psReprojInfo->poForwardTransform->GetSourceCS() == nullptr)
        return false;

    double z = 0;
    int success = true;
    auto poSourceCRS = psReprojInfo->poForwardTransform->GetSourceCS();
    if (poSourceCRS->IsGeographic() &&
        std::fabs(poSourceCRS->GetAngularUnits() -
                  CPLAtof(SRS_UA_DEGREE_CONV)) < 1e-9)
    {
        // Optimization to avoid creating a OGRCoordinateTransformation
        OGRAxisOrientation eSourceFirstAxisOrient = OAO_Other;
        poSourceCRS->GetAxis(nullptr, 0, &eSourceFirstAxisOrient);
        const auto &mapping = poSourceCRS->GetDataAxisToSRSAxisMapping();
        if ((mapping[0] == 2 && eSourceFirstAxisOrient == OAO_East) ||
            (mapping[0] == 1 && eSourceFirstAxisOrient != OAO_East))
        {
            std::swap(*pdfX, *pdfY);
        }
    }
    else
    {
        auto poLongLat =
            std::unique_ptr<OGRSpatialReference>(poSourceCRS->CloneGeogCS());
        if (poLongLat == nullptr)
            return false;
        poLongLat->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

        const bool bCurrentCheckWithInvertProj =
            GetCurrentCheckWithInvertPROJ();
        if (!bCurrentCheckWithInvertProj)
            CPLSetThreadLocalConfigOption("CHECK_WITH_INVERT_PROJ", "YES");
        auto poCT = std::unique_ptr<OGRCoordinateTransformation>(
            OGRCreateCoordinateTransformation(poLongLat.get(), poSourceCRS));
        if (!bCurrentCheckWithInvertProj)
            CPLSetThreadLocalConfigOption("CHECK_WITH_INVERT_PROJ", nullptr);
        if (poCT == nullptr)
            return false;

        poCT->SetEmitErrors(false);
        if (!poCT->Transform(1, pdfX, pdfY))
            return false;

        if (!psInfo->pReproject(psInfo->pReprojectArg, false, 1, pdfX, pdfY, &z,
                                &success) ||
            !success)
        {
            return false;
        }
    }

    double *padfGeoTransform = psInfo->sDstParams.adfInvGeoTransform;
    void *pTransformArg = psInfo->sDstParams.pTransformArg;
    GDALTransformerFunc pTransformer = psInfo->sDstParams.pTransformer;
    if (pTransformArg != nullptr)
    {
        if (!pTransformer(pTransformArg, TRUE, 1, pdfX, pdfY, &z, &success) ||
            !success)
        {
            return false;
        }
    }
    else
    {
        const double dfNewX = padfGeoTransform[0] +
                              pdfX[0] * padfGeoTransform[1] +
                              pdfY[0] * padfGeoTransform[2];
        const double dfNewY = padfGeoTransform[3] +
                              pdfX[0] * padfGeoTransform[4] +
                              pdfY[0] * padfGeoTransform[5];

        pdfX[0] = dfNewX;
        pdfY[0] = dfNewY;
    }

    return true;
}

/************************************************************************/
/*                 GDALSerializeGenImgProjTransformer()                 */
/************************************************************************/

static CPLXMLNode *GDALSerializeGenImgProjTransformer(void *pTransformArg)

{
    GDALGenImgProjTransformInfo *psInfo =
        static_cast<GDALGenImgProjTransformInfo *>(pTransformArg);

    CPLXMLNode *psTree =
        CPLCreateXMLNode(nullptr, CXT_Element, "GenImgProjTransformer");

    const auto SerializePart =
        [psTree](const char *pszPrefix, const GDALGenImgProjTransformPart &part)
    {
        char szWork[200] = {};

        /* ------------------------------------------------------------- */
        /*      Handle transformation.                                   */
        /* ------------------------------------------------------------- */
        if (part.pTransformArg != nullptr)
        {
            CPLXMLNode *psTransformer =
                GDALSerializeTransformer(part.pTransformer, part.pTransformArg);
            if (psTransformer != nullptr)
            {
                CPLXMLNode *psTransformerContainer = CPLCreateXMLNode(
                    psTree, CXT_Element,
                    CPLSPrintf("%s%s", pszPrefix, psTransformer->pszValue));

                CPLAddXMLChild(psTransformerContainer, psTransformer);
            }
        }

        /* ------------------------------------------------------------- */
        /*      Handle geotransforms.                                    */
        /* ------------------------------------------------------------- */
        else
        {
            CPLsnprintf(szWork, sizeof(szWork),
                        "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",
                        part.adfGeoTransform[0], part.adfGeoTransform[1],
                        part.adfGeoTransform[2], part.adfGeoTransform[3],
                        part.adfGeoTransform[4], part.adfGeoTransform[5]);
            CPLCreateXMLElementAndValue(
                psTree, CPLSPrintf("%sGeoTransform", pszPrefix), szWork);

            CPLsnprintf(szWork, sizeof(szWork),
                        "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g",
                        part.adfInvGeoTransform[0], part.adfInvGeoTransform[1],
                        part.adfInvGeoTransform[2], part.adfInvGeoTransform[3],
                        part.adfInvGeoTransform[4], part.adfInvGeoTransform[5]);
            CPLCreateXMLElementAndValue(
                psTree, CPLSPrintf("%sInvGeoTransform", pszPrefix), szWork);
        }
    };

    SerializePart("Src", psInfo->sSrcParams);
    SerializePart("Dst", psInfo->sDstParams);

    /* -------------------------------------------------------------------- */
    /*      Do we have a reprojection transformer?                          */
    /* -------------------------------------------------------------------- */
    if (psInfo->pReprojectArg != nullptr)
    {

        CPLXMLNode *psTransformerContainer =
            CPLCreateXMLNode(psTree, CXT_Element, "ReprojectTransformer");

        CPLXMLNode *psTransformer =
            GDALSerializeTransformer(psInfo->pReproject, psInfo->pReprojectArg);
        if (psTransformer != nullptr)
            CPLAddXMLChild(psTransformerContainer, psTransformer);
    }

    return psTree;
}

/************************************************************************/
/*                    GDALDeserializeGeoTransform()                     */
/************************************************************************/

static void GDALDeserializeGeoTransform(const char *pszGT,
                                        double adfGeoTransform[6])
{
    CPLsscanf(pszGT, "%lf,%lf,%lf,%lf,%lf,%lf", adfGeoTransform + 0,
              adfGeoTransform + 1, adfGeoTransform + 2, adfGeoTransform + 3,
              adfGeoTransform + 4, adfGeoTransform + 5);
}

/************************************************************************/
/*                GDALDeserializeGenImgProjTransformer()                */
/************************************************************************/

void *GDALDeserializeGenImgProjTransformer(CPLXMLNode *psTree)

{
    /* -------------------------------------------------------------------- */
    /*      Initialize the transform info.                                  */
    /* -------------------------------------------------------------------- */
    GDALGenImgProjTransformInfo *psInfo =
        GDALCreateGenImgProjTransformerInternal();

    const auto DeserializePart =
        [psTree](const char *pszPrefix, GDALGenImgProjTransformPart &part)
    {
        /* ----------------------------------------------------------------- */
        /*      Geotransform                                                 */
        /* ----------------------------------------------------------------- */
        if (const auto psGTNode =
                CPLGetXMLNode(psTree, CPLSPrintf("%sGeoTransform", pszPrefix)))
        {
            GDALDeserializeGeoTransform(CPLGetXMLValue(psGTNode, "", ""),
                                        part.adfGeoTransform);

            if (const auto psInvGTNode = CPLGetXMLNode(
                    psTree, CPLSPrintf("%sInvGeoTransform", pszPrefix)))
            {
                GDALDeserializeGeoTransform(CPLGetXMLValue(psInvGTNode, "", ""),
                                            part.adfInvGeoTransform);
            }
            else
            {
                if (!GDALInvGeoTransform(part.adfGeoTransform,
                                         part.adfInvGeoTransform))
                {
                    CPLError(CE_Failure, CPLE_AppDefined,
                             "Cannot invert geotransform");
                }
            }
        }

        /* ---------------------------------------------------------------- */
        /*      Transform                                                   */
        /* ---------------------------------------------------------------- */
        else
        {
            for (CPLXMLNode *psIter = psTree->psChild; psIter != nullptr;
                 psIter = psIter->psNext)
            {
                if (psIter->eType == CXT_Element &&
                    STARTS_WITH_CI(psIter->pszValue, pszPrefix))
                {
                    GDALDeserializeTransformer(psIter->psChild,
                                               &part.pTransformer,
                                               &part.pTransformArg);
                    break;
                }
            }
        }
    };

    DeserializePart("Src", psInfo->sSrcParams);
    DeserializePart("Dst", psInfo->sDstParams);

    /* -------------------------------------------------------------------- */
    /*      Reproject transformer                                           */
    /* -------------------------------------------------------------------- */
    CPLXMLNode *psSubtree = CPLGetXMLNode(psTree, "ReprojectTransformer");
    if (psSubtree != nullptr && psSubtree->psChild != nullptr)
    {
        GDALDeserializeTransformer(psSubtree->psChild, &psInfo->pReproject,
                                   &psInfo->pReprojectArg);
    }

    return psInfo;
}

/************************************************************************/
/*                 GDALCreateReprojectionTransformer()                  */
/************************************************************************/

/**
 * Create reprojection transformer.
 *
 * Creates a callback data structure suitable for use with
 * GDALReprojectionTransformation() to represent a transformation from
 * one geographic or projected coordinate system to another.  On input
 * the coordinate systems are described in OpenGIS WKT format.
 *
 * Internally the OGRCoordinateTransformation object is used to implement
 * the reprojection.
 *
 * @param pszSrcWKT the coordinate system for the source coordinate system.
 * @param pszDstWKT the coordinate system for the destination coordinate
 * system.
 *
 * @return Handle for use with GDALReprojectionTransform(), or NULL if the
 * system fails to initialize the reprojection.
 **/

void *GDALCreateReprojectionTransformer(const char *pszSrcWKT,
                                        const char *pszDstWKT)

{
    /* -------------------------------------------------------------------- */
    /*      Ingest the SRS definitions.                                     */
    /* -------------------------------------------------------------------- */
    OGRSpatialReference oSrcSRS;
    oSrcSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (oSrcSRS.importFromWkt(pszSrcWKT) != OGRERR_NONE)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Failed to import coordinate system `%s'.", pszSrcWKT);
        return nullptr;
    }

    OGRSpatialReference oDstSRS;
    oDstSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (oDstSRS.importFromWkt(pszDstWKT) != OGRERR_NONE)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Failed to import coordinate system `%s'.", pszSrcWKT);
        return nullptr;
    }

    return GDALCreateReprojectionTransformerEx(
        OGRSpatialReference::ToHandle(&oSrcSRS),
        OGRSpatialReference::ToHandle(&oDstSRS), nullptr);
}

/************************************************************************/
/*                 GDALCreateReprojectionTransformerEx()                */
/************************************************************************/

/**
 * Create reprojection transformer.
 *
 * Creates a callback data structure suitable for use with
 * GDALReprojectionTransformation() to represent a transformation from
 * one geographic or projected coordinate system to another.
 *
 * Internally the OGRCoordinateTransformation object is used to implement
 * the reprojection.
 *
 * @param hSrcSRS the coordinate system for the source coordinate system.
 * @param hDstSRS the coordinate system for the destination coordinate
 * system.
 * @param papszOptions NULL-terminated list of options, or NULL. Currently
 * supported options are:
 * <ul>
 * <li>AREA_OF_INTEREST=west_long,south_lat,east_long,north_lat: Values in
 * degrees. longitudes in [-180,180], latitudes in [-90,90].</li>
 * <li>COORDINATE_OPERATION=string: PROJ or WKT string representing a
 * coordinate operation, overriding the default computed transformation.</li>
 * <li>COORDINATE_EPOCH=decimal_year: Coordinate epoch, expressed as a
 * decimal year. Useful for time-dependent coordinate operations.</li>
 * <li> SRC_COORDINATE_EPOCH: (GDAL &gt;= 3.4) Coordinate epoch of source CRS,
 * expressed as a decimal year. Useful for time-dependent coordinate
 *operations.</li>
 * <li> DST_COORDINATE_EPOCH: (GDAL &gt;= 3.4) Coordinate epoch
 *of target CRS, expressed as a decimal year. Useful for time-dependent
 *coordinate operations.</li>
 * <li> ALLOW_BALLPARK=YES/NO: (GDAL &gt;= 3.11) Whether ballpark coordinate
 * operations are allowed. Defaults to YES.</li>
 * <li> ONLY_BEST=YES/NO/AUTO: (GDAL &gt;= 3.11) By default (at least in the
 * PROJ 9.x series), PROJ may use coordinate
 * operations that are not the "best" if resources (typically grids) needed
 * to use them are missing. It will then fallback to other coordinate operations
 * that have a lesser accuracy, for example using Helmert transformations,
 * or in the absence of such operations, to ones with potential very rought
 * accuracy, using "ballpark" transformations
 * (see https://proj.org/glossary.html).
 * When calling this method with YES, PROJ will only consider the
 * "best" operation, and error out (at Transform() time) if they cannot be
 * used.
 * This method may be used together with ALLOW_BALLPARK=NO to
 * only allow best operations that have a known accuracy.
 * Note that this method has no effect on PROJ versions before 9.2.
 * The default value for this option can be also set with the
 * PROJ_ONLY_BEST_DEFAULT environment variable, or with the "only_best_default"
 * setting of proj.ini. Calling SetOnlyBest() overrides such default value.</li>
 * </ul>
 *
 * @return Handle for use with GDALReprojectionTransform(), or NULL if the
 * system fails to initialize the reprojection.
 *
 * @since GDAL 3.0
 **/

void *GDALCreateReprojectionTransformerEx(OGRSpatialReferenceH hSrcSRS,
                                          OGRSpatialReferenceH hDstSRS,
                                          const char *const *papszOptions)
{
    OGRSpatialReference *poSrcSRS = OGRSpatialReference::FromHandle(hSrcSRS);
    OGRSpatialReference *poDstSRS = OGRSpatialReference::FromHandle(hDstSRS);

    /* -------------------------------------------------------------------- */
    /*      Build the forward coordinate transformation.                    */
    /* -------------------------------------------------------------------- */
    double dfWestLongitudeDeg = 0.0;
    double dfSouthLatitudeDeg = 0.0;
    double dfEastLongitudeDeg = 0.0;
    double dfNorthLatitudeDeg = 0.0;
    const char *pszBBOX = CSLFetchNameValue(papszOptions, "AREA_OF_INTEREST");
    if (pszBBOX)
    {
        char **papszTokens = CSLTokenizeString2(pszBBOX, ",", 0);
        if (CSLCount(papszTokens) == 4)
        {
            dfWestLongitudeDeg = CPLAtof(papszTokens[0]);
            dfSouthLatitudeDeg = CPLAtof(papszTokens[1]);
            dfEastLongitudeDeg = CPLAtof(papszTokens[2]);
            dfNorthLatitudeDeg = CPLAtof(papszTokens[3]);
        }
        CSLDestroy(papszTokens);
    }
    const char *pszCO = CSLFetchNameValue(papszOptions, "COORDINATE_OPERATION");

    OGRCoordinateTransformationOptions optionsFwd;
    if (!(dfWestLongitudeDeg == 0.0 && dfSouthLatitudeDeg == 0.0 &&
          dfEastLongitudeDeg == 0.0 && dfNorthLatitudeDeg == 0.0))
    {
        optionsFwd.SetAreaOfInterest(dfWestLongitudeDeg, dfSouthLatitudeDeg,
                                     dfEastLongitudeDeg, dfNorthLatitudeDeg);
    }
    if (pszCO)
    {
        optionsFwd.SetCoordinateOperation(pszCO, false);
    }

    const char *pszCENTER_LONG = CSLFetchNameValue(papszOptions, "CENTER_LONG");
    if (pszCENTER_LONG)
    {
        optionsFwd.SetSourceCenterLong(CPLAtof(pszCENTER_LONG));
    }

    optionsFwd.SetBallparkAllowed(CPLTestBool(
        CSLFetchNameValueDef(papszOptions, "ALLOW_BALLPARK", "YES")));

    const char *pszOnlyBest =
        CSLFetchNameValueDef(papszOptions, "ONLY_BEST", "AUTO");
    if (!EQUAL(pszOnlyBest, "AUTO"))
    {
        optionsFwd.SetOnlyBest(CPLTestBool(pszOnlyBest));
    }

    OGRCoordinateTransformation *poForwardTransform =
        OGRCreateCoordinateTransformation(poSrcSRS, poDstSRS, optionsFwd);

    if (poForwardTransform == nullptr)
        // OGRCreateCoordinateTransformation() will report errors on its own.
        return nullptr;

    poForwardTransform->SetEmitErrors(false);

    /* -------------------------------------------------------------------- */
    /*      Create a structure to hold the transform info, and also         */
    /*      build reverse transform.  We assume that if the forward         */
    /*      transform can be created, then so can the reverse one.          */
    /* -------------------------------------------------------------------- */
    GDALReprojectionTransformInfo *psInfo = new GDALReprojectionTransformInfo();

    psInfo->papszOptions = CSLDuplicate(papszOptions);
    psInfo->poForwardTransform = poForwardTransform;
    psInfo->dfTime = CPLAtof(CSLFetchNameValueDef(
        papszOptions, "COORDINATE_EPOCH",
        CSLFetchNameValueDef(
            papszOptions, "DST_COORDINATE_EPOCH",
            CSLFetchNameValueDef(papszOptions, "SRC_COORDINATE_EPOCH", "0"))));
    psInfo->poReverseTransform = poForwardTransform->GetInverse();

    if (psInfo->poReverseTransform)
        psInfo->poReverseTransform->SetEmitErrors(false);

    memcpy(psInfo->sTI.abySignature, GDAL_GTI2_SIGNATURE,
           strlen(GDAL_GTI2_SIGNATURE));
    psInfo->sTI.pszClassName = GDAL_REPROJECTION_TRANSFORMER_CLASS_NAME;
    psInfo->sTI.pfnTransform = GDALReprojectionTransform;
    psInfo->sTI.pfnCleanup = GDALDestroyReprojectionTransformer;
    psInfo->sTI.pfnSerialize = GDALSerializeReprojectionTransformer;

    return psInfo;
}

/************************************************************************/
/*                 GDALDestroyReprojectionTransformer()                 */
/************************************************************************/

/**
 * Destroy reprojection transformation.
 *
 * @param pTransformArg the transformation handle returned by
 * GDALCreateReprojectionTransformer().
 */

void GDALDestroyReprojectionTransformer(void *pTransformArg)

{
    if (pTransformArg == nullptr)
        return;

    GDALReprojectionTransformInfo *psInfo =
        static_cast<GDALReprojectionTransformInfo *>(pTransformArg);

    if (psInfo->poForwardTransform)
        OGRCoordinateTransformation::DestroyCT(psInfo->poForwardTransform);

    if (psInfo->poReverseTransform)
        OGRCoordinateTransformation::DestroyCT(psInfo->poReverseTransform);

    CSLDestroy(psInfo->papszOptions);

    delete psInfo;
}

/************************************************************************/
/*                     GDALReprojectionTransform()                      */
/************************************************************************/

/**
 * Perform reprojection transformation.
 *
 * Actually performs the reprojection transformation described in
 * GDALCreateReprojectionTransformer().  This function matches the
 * GDALTransformerFunc() signature.  Details of the arguments are described
 * there.
 */

int GDALReprojectionTransform(void *pTransformArg, int bDstToSrc,
                              int nPointCount, double *padfX, double *padfY,
                              double *padfZ, int *panSuccess)

{
    GDALReprojectionTransformInfo *psInfo =
        static_cast<GDALReprojectionTransformInfo *>(pTransformArg);
    int bSuccess;

    std::vector<double> adfTime;
    double *padfT = nullptr;
    if (psInfo->dfTime != 0.0 && nPointCount > 0)
    {
        adfTime.resize(nPointCount, psInfo->dfTime);
        padfT = &adfTime[0];
    }

    if (bDstToSrc)
    {
        if (psInfo->poReverseTransform == nullptr)
        {
            CPLError(
                CE_Failure, CPLE_AppDefined,
                "Inverse coordinate transformation cannot be instantiated");
            if (panSuccess)
            {
                for (int i = 0; i < nPointCount; i++)
                    panSuccess[i] = FALSE;
            }
            bSuccess = false;
        }
        else
        {
            bSuccess = psInfo->poReverseTransform->Transform(
                nPointCount, padfX, padfY, padfZ, padfT, panSuccess);
        }
    }
    else
        bSuccess = psInfo->poForwardTransform->Transform(
            nPointCount, padfX, padfY, padfZ, padfT, panSuccess);

    return bSuccess;
}

/************************************************************************/
/*                GDALSerializeReprojectionTransformer()                */
/************************************************************************/

static CPLXMLNode *GDALSerializeReprojectionTransformer(void *pTransformArg)

{
    CPLXMLNode *psTree;
    GDALReprojectionTransformInfo *psInfo =
        static_cast<GDALReprojectionTransformInfo *>(pTransformArg);

    psTree = CPLCreateXMLNode(nullptr, CXT_Element, "ReprojectionTransformer");

    /* -------------------------------------------------------------------- */
    /*      Handle SourceCS.                                                */
    /* -------------------------------------------------------------------- */
    const auto ExportToWkt = [](const OGRSpatialReference *poSRS)
    {
        // Try first in WKT1 for backward compat
        {
            char *pszWKT = nullptr;
            const char *const apszOptions[] = {"FORMAT=WKT1", nullptr};
            CPLErrorHandlerPusher oHandler(CPLQuietErrorHandler);
            CPLErrorStateBackuper oBackuper;
            if (poSRS->exportToWkt(&pszWKT, apszOptions) == OGRERR_NONE)
            {
                std::string osRet(pszWKT);
                CPLFree(pszWKT);
                return osRet;
            }
            CPLFree(pszWKT);
        }

        char *pszWKT = nullptr;
        const char *const apszOptions[] = {"FORMAT=WKT2_2019", nullptr};
        if (poSRS->exportToWkt(&pszWKT, apszOptions) == OGRERR_NONE)
        {
            std::string osRet(pszWKT);
            CPLFree(pszWKT);
            return osRet;
        }
        CPLFree(pszWKT);
        return std::string();
    };

    auto poSRS = psInfo->poForwardTransform->GetSourceCS();
    if (poSRS)
    {
        const auto osWKT = ExportToWkt(poSRS);
        CPLCreateXMLElementAndValue(psTree, "SourceSRS", osWKT.c_str());
    }

    /* -------------------------------------------------------------------- */
    /*      Handle DestinationCS.                                           */
    /* -------------------------------------------------------------------- */
    poSRS = psInfo->poForwardTransform->GetTargetCS();
    if (poSRS)
    {
        const auto osWKT = ExportToWkt(poSRS);
        CPLCreateXMLElementAndValue(psTree, "TargetSRS", osWKT.c_str());
    }

    /* -------------------------------------------------------------------- */
    /*      Serialize options.                                              */
    /* -------------------------------------------------------------------- */
    if (psInfo->papszOptions)
    {
        CPLXMLNode *psOptions =
            CPLCreateXMLNode(psTree, CXT_Element, "Options");
        for (auto iter = psInfo->papszOptions; *iter != nullptr; ++iter)
        {
            char *pszKey = nullptr;
            const char *pszValue = CPLParseNameValue(*iter, &pszKey);
            if (pszKey && pszValue)
            {
                auto elt =
                    CPLCreateXMLElementAndValue(psOptions, "Option", pszValue);
                CPLAddXMLAttributeAndValue(elt, "key", pszKey);
            }
            CPLFree(pszKey);
        }
    }

    return psTree;
}

/************************************************************************/
/*               GDALDeserializeReprojectionTransformer()               */
/************************************************************************/

static void *GDALDeserializeReprojectionTransformer(CPLXMLNode *psTree)

{
    const char *pszSourceSRS = CPLGetXMLValue(psTree, "SourceSRS", nullptr);
    const char *pszTargetSRS = CPLGetXMLValue(psTree, "TargetSRS", nullptr);
    char *pszSourceWKT = nullptr, *pszTargetWKT = nullptr;
    void *pResult = nullptr;

    OGRSpatialReference oSrcSRS;
    OGRSpatialReference oDstSRS;

    oSrcSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    oDstSRS.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (pszSourceSRS != nullptr)
    {
        oSrcSRS.SetFromUserInput(pszSourceSRS);
    }

    if (pszTargetSRS != nullptr)
    {
        oDstSRS.SetFromUserInput(pszTargetSRS);
    }

    CPLStringList aosList;
    const CPLXMLNode *psOptions = CPLGetXMLNode(psTree, "Options");
    if (psOptions)
    {
        for (auto iter = psOptions->psChild; iter; iter = iter->psNext)
        {
            if (iter->eType == CXT_Element &&
                strcmp(iter->pszValue, "Option") == 0)
            {
                const char *pszKey = CPLGetXMLValue(iter, "key", nullptr);
                const char *pszValue = CPLGetXMLValue(iter, nullptr, nullptr);
                if (pszKey && pszValue)
                {
                    aosList.SetNameValue(pszKey, pszValue);
                }
            }
        }
    }

    pResult = GDALCreateReprojectionTransformerEx(
        !oSrcSRS.IsEmpty() ? OGRSpatialReference::ToHandle(&oSrcSRS) : nullptr,
        !oDstSRS.IsEmpty() ? OGRSpatialReference::ToHandle(&oDstSRS) : nullptr,
        aosList.List());

    CPLFree(pszSourceWKT);
    CPLFree(pszTargetWKT);

    return pResult;
}

/************************************************************************/
/* ==================================================================== */
/*      Approximate transformer.                                        */
/* ==================================================================== */
/************************************************************************/

/************************************************************************/
/*                  GDALCreateSimilarApproxTransformer()                */
/************************************************************************/

static void *GDALCreateSimilarApproxTransformer(void *hTransformArg,
                                                double dfSrcRatioX,
                                                double dfSrcRatioY)
{
    VALIDATE_POINTER1(hTransformArg, "GDALCreateSimilarApproxTransformer",
                      nullptr);

    GDALApproxTransformInfo *psInfo =
        static_cast<GDALApproxTransformInfo *>(hTransformArg);

    void *pBaseCBData = GDALCreateSimilarTransformer(psInfo->pBaseCBData,
                                                     dfSrcRatioX, dfSrcRatioY);
    if (pBaseCBData == nullptr)
    {
        return nullptr;
    }

    GDALApproxTransformInfo *psClonedInfo =
        static_cast<GDALApproxTransformInfo *>(GDALCreateApproxTransformer2(
            psInfo->pfnBaseTransformer, pBaseCBData, psInfo->dfMaxErrorForward,
            psInfo->dfMaxErrorReverse));
    psClonedInfo->bOwnSubtransformer = TRUE;

    return psClonedInfo;
}

/************************************************************************/
/*                   GDALSerializeApproxTransformer()                   */
/************************************************************************/

static CPLXMLNode *GDALSerializeApproxTransformer(void *pTransformArg)

{
    CPLXMLNode *psTree;
    GDALApproxTransformInfo *psInfo =
        static_cast<GDALApproxTransformInfo *>(pTransformArg);

    psTree = CPLCreateXMLNode(nullptr, CXT_Element, "ApproxTransformer");

    /* -------------------------------------------------------------------- */
    /*      Attach max error.                                               */
    /* -------------------------------------------------------------------- */
    if (psInfo->dfMaxErrorForward == psInfo->dfMaxErrorReverse)
    {
        CPLCreateXMLElementAndValue(
            psTree, "MaxError",
            CPLString().Printf("%g", psInfo->dfMaxErrorForward));
    }
    else
    {
        CPLCreateXMLElementAndValue(
            psTree, "MaxErrorForward",
            CPLString().Printf("%g", psInfo->dfMaxErrorForward));
        CPLCreateXMLElementAndValue(
            psTree, "MaxErrorReverse",
            CPLString().Printf("%g", psInfo->dfMaxErrorReverse));
    }

    /* -------------------------------------------------------------------- */
    /*      Capture underlying transformer.                                 */
    /* -------------------------------------------------------------------- */
    CPLXMLNode *psTransformerContainer =
        CPLCreateXMLNode(psTree, CXT_Element, "BaseTransformer");

    CPLXMLNode *psTransformer = GDALSerializeTransformer(
        psInfo->pfnBaseTransformer, psInfo->pBaseCBData);
    if (psTransformer != nullptr)
        CPLAddXMLChild(psTransformerContainer, psTransformer);

    return psTree;
}

/************************************************************************/
/*                    GDALCreateApproxTransformer()                     */
/************************************************************************/

/**
 * Create an approximating transformer.
 *
 * This function creates a context for an approximated transformer.  Basically
 * a high precision transformer is supplied as input and internally linear
 * approximations are computed to generate results to within a defined
 * precision.
 *
 * The approximation is actually done at the point where GDALApproxTransform()
 * calls are made, and depend on the assumption that they are roughly linear.
 * The first and last point passed in must be the extreme values and the
 * intermediate values should describe a curve between the end points.  The
 * approximator transforms and centers using the approximate transformer, and
 * then compares the true middle transformed value to a linear approximation
 * based on the end points.  If the error is within the supplied threshold then
 * the end points are used to linearly approximate all the values otherwise the
 * input points are split into two smaller sets, and the function is recursively
 * called until a sufficiently small set of points is found that the linear
 * approximation is OK, or that all the points are exactly computed.
 *
 * This function is very suitable for approximating transformation results
 * from output pixel/line space to input coordinates for warpers that operate
 * on one input scanline at a time.  Care should be taken using it in other
 * circumstances as little internal validation is done in order to keep things
 * fast.
 *
 * @param pfnBaseTransformer the high precision transformer which should be
 * approximated.
 * @param pBaseTransformArg the callback argument for the high precision
 * transformer.
 * @param dfMaxError the maximum cartesian error in the "output" space that
 * is to be accepted in the linear approximation, evaluated as a Manhattan
 * distance.
 *
 * @return callback pointer suitable for use with GDALApproxTransform().  It
 * should be deallocated with GDALDestroyApproxTransformer().
 */

void *GDALCreateApproxTransformer(GDALTransformerFunc pfnBaseTransformer,
                                  void *pBaseTransformArg, double dfMaxError)

{
    return GDALCreateApproxTransformer2(pfnBaseTransformer, pBaseTransformArg,
                                        dfMaxError, dfMaxError);
}

static void *
GDALCreateApproxTransformer2(GDALTransformerFunc pfnBaseTransformer,
                             void *pBaseTransformArg, double dfMaxErrorForward,
                             double dfMaxErrorReverse)

{
    GDALApproxTransformInfo *psATInfo = new GDALApproxTransformInfo;
    psATInfo->pfnBaseTransformer = pfnBaseTransformer;
    psATInfo->pBaseCBData = pBaseTransformArg;
    psATInfo->dfMaxErrorForward = dfMaxErrorForward;
    psATInfo->dfMaxErrorReverse = dfMaxErrorReverse;
    psATInfo->bOwnSubtransformer = FALSE;

    memcpy(psATInfo->sTI.abySignature, GDAL_GTI2_SIGNATURE,
           strlen(GDAL_GTI2_SIGNATURE));
    psATInfo->sTI.pszClassName = GDAL_APPROX_TRANSFORMER_CLASS_NAME;
    psATInfo->sTI.pfnTransform = GDALApproxTransform;
    psATInfo->sTI.pfnCleanup = GDALDestroyApproxTransformer;
    psATInfo->sTI.pfnSerialize = GDALSerializeApproxTransformer;
    psATInfo->sTI.pfnCreateSimilar = GDALCreateSimilarApproxTransformer;

    return psATInfo;
}

/************************************************************************/
/*              GDALApproxTransformerOwnsSubtransformer()               */
/************************************************************************/

/** Set bOwnSubtransformer flag */
void GDALApproxTransformerOwnsSubtransformer(void *pCBData, int bOwnFlag)

{
    GDALApproxTransformInfo *psATInfo =
        static_cast<GDALApproxTransformInfo *>(pCBData);

    psATInfo->bOwnSubtransformer = bOwnFlag;
}

/************************************************************************/
/*                    GDALDestroyApproxTransformer()                    */
/************************************************************************/

/**
 * Cleanup approximate transformer.
 *
 * Deallocates the resources allocated by GDALCreateApproxTransformer().
 *
 * @param pCBData callback data originally returned by
 * GDALCreateApproxTransformer().
 */

void GDALDestroyApproxTransformer(void *pCBData)

{
    if (pCBData == nullptr)
        return;

    GDALApproxTransformInfo *psATInfo =
        static_cast<GDALApproxTransformInfo *>(pCBData);

    if (psATInfo->bOwnSubtransformer)
        GDALDestroyTransformer(psATInfo->pBaseCBData);

    delete psATInfo;
}

/************************************************************************/
/*                  GDALRefreshApproxTransformer()                      */
/************************************************************************/

void GDALRefreshApproxTransformer(void *hTransformArg)
{
    GDALApproxTransformInfo *psInfo =
        static_cast<GDALApproxTransformInfo *>(hTransformArg);

    if (GDALIsTransformer(psInfo->pBaseCBData,
                          GDAL_GEN_IMG_TRANSFORMER_CLASS_NAME))
    {
        GDALRefreshGenImgProjTransformer(psInfo->pBaseCBData);
    }
}

/************************************************************************/
/*                      GDALApproxTransformInternal()                   */
/************************************************************************/

static int GDALApproxTransformInternal(void *pCBData, int bDstToSrc,
                                       int nPoints, double *x, double *y,
                                       double *z, int *panSuccess,
                                       // SME = Start, Middle, End.
                                       const double xSMETransformed[3],
                                       const double ySMETransformed[3],
                                       const double zSMETransformed[3])
{
    GDALApproxTransformInfo *psATInfo =
        static_cast<GDALApproxTransformInfo *>(pCBData);
    const int nMiddle = (nPoints - 1) / 2;

#ifdef notdef_sanify_check
    {
        double x2[3] = {x[0], x[nMiddle], x[nPoints - 1]};
        double y2[3] = {y[0], y[nMiddle], y[nPoints - 1]};
        double z2[3] = {z[0], z[nMiddle], z[nPoints - 1]};
        int anSuccess2[3] = {};

        const int bSuccess = psATInfo->pfnBaseTransformer(
            psATInfo->pBaseCBData, bDstToSrc, 3, x2, y2, z2, anSuccess2);
        CPLAssert(bSuccess);
        CPLAssert(anSuccess2[0]);
        CPLAssert(anSuccess2[1]);
        CPLAssert(anSuccess2[2]);
        CPLAssert(x2[0] == xSMETransformed[0]);
        CPLAssert(y2[0] == ySMETransformed[0]);
        CPLAssert(z2[0] == zSMETransformed[0]);
        CPLAssert(x2[1] == xSMETransformed[1]);
        CPLAssert(y2[1] == ySMETransformed[1]);
        CPLAssert(z2[1] == zSMETransformed[1]);
        CPLAssert(x2[2] == xSMETransformed[2]);
        CPLAssert(y2[2] == ySMETransformed[2]);
        CPLAssert(z2[2] == zSMETransformed[2]);
    }
#endif

#ifdef DEBUG_APPROX_TRANSFORMER
    fprintf(stderr, "start (%.3f,%.3f) -> (%.3f,%.3f)\n", /*ok*/
            x[0], y[0], xSMETransformed[0], ySMETransformed[0]);
    fprintf(stderr, "middle (%.3f,%.3f) -> (%.3f,%.3f)\n", /*ok*/
            x[nMiddle], y[nMiddle], xSMETransformed[1], ySMETransformed[1]);
    fprintf(stderr, "end (%.3f,%.3f) -> (%.3f,%.3f)\n", /*ok*/
            x[nPoints - 1], y[nPoints - 1], xSMETransformed[2],
            ySMETransformed[2]);
#endif

    /* -------------------------------------------------------------------- */
    /*      Is the error at the middle acceptable relative to an            */
    /*      interpolation of the middle position?                           */
    /* -------------------------------------------------------------------- */
    const double dfDeltaX =
        (xSMETransformed[2] - xSMETransformed[0]) / (x[nPoints - 1] - x[0]);
    const double dfDeltaY =
        (ySMETransformed[2] - ySMETransformed[0]) / (x[nPoints - 1] - x[0]);
    const double dfDeltaZ =
        (zSMETransformed[2] - zSMETransformed[0]) / (x[nPoints - 1] - x[0]);

    const double dfError =
        fabs((xSMETransformed[0] + dfDeltaX * (x[nMiddle] - x[0])) -
             xSMETransformed[1]) +
        fabs((ySMETransformed[0] + dfDeltaY * (x[nMiddle] - x[0])) -
             ySMETransformed[1]);

    const double dfMaxError =
        (bDstToSrc) ? psATInfo->dfMaxErrorReverse : psATInfo->dfMaxErrorForward;
    if (dfError > dfMaxError)
    {
#if DEBUG_VERBOSE
        CPLDebug("GDAL",
                 "ApproxTransformer - "
                 "error %g over threshold %g, subdivide %d points.",
                 dfError, dfMaxError, nPoints);
#endif

        double xMiddle[3] = {x[(nMiddle - 1) / 2], x[nMiddle - 1],
                             x[nMiddle + (nPoints - nMiddle - 1) / 2]};
        double yMiddle[3] = {y[(nMiddle - 1) / 2], y[nMiddle - 1],
                             y[nMiddle + (nPoints - nMiddle - 1) / 2]};
        double zMiddle[3] = {z[(nMiddle - 1) / 2], z[nMiddle - 1],
                             z[nMiddle + (nPoints - nMiddle - 1) / 2]};

        const bool bUseBaseTransformForHalf1 =
            nMiddle <= 5 || y[0] != y[nMiddle - 1] ||
            y[0] != y[(nMiddle - 1) / 2] || x[0] == x[nMiddle - 1] ||
            x[0] == x[(nMiddle - 1) / 2];
        const bool bUseBaseTransformForHalf2 =
            nPoints - nMiddle <= 5 || y[nMiddle] != y[nPoints - 1] ||
            y[nMiddle] != y[nMiddle + (nPoints - nMiddle - 1) / 2] ||
            x[nMiddle] == x[nPoints - 1] ||
            x[nMiddle] == x[nMiddle + (nPoints - nMiddle - 1) / 2];

        int anSuccess2[3] = {};
        int bSuccess = FALSE;
        if (!bUseBaseTransformForHalf1 && !bUseBaseTransformForHalf2)
            bSuccess = psATInfo->pfnBaseTransformer(
                psATInfo->pBaseCBData, bDstToSrc, 3, xMiddle, yMiddle, zMiddle,
                anSuccess2);
        else if (!bUseBaseTransformForHalf1)
        {
            bSuccess = psATInfo->pfnBaseTransformer(
                psATInfo->pBaseCBData, bDstToSrc, 2, xMiddle, yMiddle, zMiddle,
                anSuccess2);
            anSuccess2[2] = TRUE;
        }
        else if (!bUseBaseTransformForHalf2)
        {
            bSuccess = psATInfo->pfnBaseTransformer(
                psATInfo->pBaseCBData, bDstToSrc, 1, xMiddle + 2, yMiddle + 2,
                zMiddle + 2, anSuccess2 + 2);
            anSuccess2[0] = TRUE;
            anSuccess2[1] = TRUE;
        }

        if (!bSuccess || !anSuccess2[0] || !anSuccess2[1] || !anSuccess2[2])
        {
            bSuccess = psATInfo->pfnBaseTransformer(
                psATInfo->pBaseCBData, bDstToSrc, nMiddle - 1, x + 1, y + 1,
                z + 1, panSuccess + 1);
            bSuccess &= psATInfo->pfnBaseTransformer(
                psATInfo->pBaseCBData, bDstToSrc, nPoints - nMiddle - 2,
                x + nMiddle + 1, y + nMiddle + 1, z + nMiddle + 1,
                panSuccess + nMiddle + 1);

            x[0] = xSMETransformed[0];
            y[0] = ySMETransformed[0];
            z[0] = zSMETransformed[0];
            panSuccess[0] = TRUE;
            x[nMiddle] = xSMETransformed[1];
            y[nMiddle] = ySMETransformed[1];
            z[nMiddle] = zSMETransformed[1];
            panSuccess[nMiddle] = TRUE;
            x[nPoints - 1] = xSMETransformed[2];
            y[nPoints - 1] = ySMETransformed[2];
            z[nPoints - 1] = zSMETransformed[2];
            panSuccess[nPoints - 1] = TRUE;
            return bSuccess;
        }

        double x2[3] = {};
        double y2[3] = {};
        double z2[3] = {};
        if (!bUseBaseTransformForHalf1)
        {
            x2[0] = xSMETransformed[0];
            y2[0] = ySMETransformed[0];
            z2[0] = zSMETransformed[0];
            x2[1] = xMiddle[0];
            y2[1] = yMiddle[0];
            z2[1] = zMiddle[0];
            x2[2] = xMiddle[1];
            y2[2] = yMiddle[1];
            z2[2] = zMiddle[1];

            bSuccess = GDALApproxTransformInternal(
                psATInfo, bDstToSrc, nMiddle, x, y, z, panSuccess, x2, y2, z2);
        }
        else
        {
            bSuccess = psATInfo->pfnBaseTransformer(
                psATInfo->pBaseCBData, bDstToSrc, nMiddle - 1, x + 1, y + 1,
                z + 1, panSuccess + 1);
            x[0] = xSMETransformed[0];
            y[0] = ySMETransformed[0];
            z[0] = zSMETransformed[0];
            panSuccess[0] = TRUE;
        }

        if (!bSuccess)
            return FALSE;

        if (!bUseBaseTransformForHalf2)
        {
            x2[0] = xSMETransformed[1];
            y2[0] = ySMETransformed[1];
            z2[0] = zSMETransformed[1];
            x2[1] = xMiddle[2];
            y2[1] = yMiddle[2];
            z2[1] = zMiddle[2];
            x2[2] = xSMETransformed[2];
            y2[2] = ySMETransformed[2];
            z2[2] = zSMETransformed[2];

            bSuccess = GDALApproxTransformInternal(
                psATInfo, bDstToSrc, nPoints - nMiddle, x + nMiddle,
                y + nMiddle, z + nMiddle, panSuccess + nMiddle, x2, y2, z2);
        }
        else
        {
            bSuccess = psATInfo->pfnBaseTransformer(
                psATInfo->pBaseCBData, bDstToSrc, nPoints - nMiddle - 2,
                x + nMiddle + 1, y + nMiddle + 1, z + nMiddle + 1,
                panSuccess + nMiddle + 1);

            x[nMiddle] = xSMETransformed[1];
            y[nMiddle] = ySMETransformed[1];
            z[nMiddle] = zSMETransformed[1];
            panSuccess[nMiddle] = TRUE;
            x[nPoints - 1] = xSMETransformed[2];
            y[nPoints - 1] = ySMETransformed[2];
            z[nPoints - 1] = zSMETransformed[2];
            panSuccess[nPoints - 1] = TRUE;
        }

        if (!bSuccess)
            return FALSE;

        return TRUE;
    }

    /* -------------------------------------------------------------------- */
    /*      Error is OK since this is just used to compute output bounds    */
    /*      of newly created file for gdalwarper.  So just use affine       */
    /*      approximation of the reverse transform.  Eventually we          */
    /*      should implement iterative searching to find a result within    */
    /*      our error threshold.                                            */
    /*      NOTE: the above comment is not true: gdalwarp uses approximator */
    /*      also to compute the source pixel of each target pixel.          */
    /* -------------------------------------------------------------------- */
    for (int i = nPoints - 1; i >= 0; i--)
    {
#ifdef check_error
        double xtemp = x[i];
        double ytemp = y[i];
        double ztemp = z[i];
        double x_ori = xtemp;
        double y_ori = ytemp;
        int btemp = FALSE;
        psATInfo->pfnBaseTransformer(psATInfo->pBaseCBData, bDstToSrc, 1,
                                     &xtemp, &ytemp, &ztemp, &btemp);
#endif
        const double dfDist = (x[i] - x[0]);
        x[i] = xSMETransformed[0] + dfDeltaX * dfDist;
        y[i] = ySMETransformed[0] + dfDeltaY * dfDist;
        z[i] = zSMETransformed[0] + dfDeltaZ * dfDist;
#ifdef check_error
        const double dfError2 = fabs(x[i] - xtemp) + fabs(y[i] - ytemp);
        if (dfError2 > 4 /*10 * dfMaxError*/)
        {
            /*ok*/ printf("Error = %f on (%f, %f)\n", dfError2, x_ori, y_ori);
        }
#endif
        panSuccess[i] = TRUE;
    }

    return TRUE;
}

/************************************************************************/
/*                        GDALApproxTransform()                         */
/************************************************************************/

/**
 * Perform approximate transformation.
 *
 * Actually performs the approximate transformation described in
 * GDALCreateApproxTransformer().  This function matches the
 * GDALTransformerFunc() signature.  Details of the arguments are described
 * there.
 */

int GDALApproxTransform(void *pCBData, int bDstToSrc, int nPoints, double *x,
                        double *y, double *z, int *panSuccess)

{
    GDALApproxTransformInfo *psATInfo =
        static_cast<GDALApproxTransformInfo *>(pCBData);
    double x2[3] = {};
    double y2[3] = {};
    double z2[3] = {};
    int anSuccess2[3] = {};
    int bSuccess;

    const int nMiddle = (nPoints - 1) / 2;

    /* -------------------------------------------------------------------- */
    /*      Bail if our preconditions are not met, or if error is not       */
    /*      acceptable.                                                     */
    /* -------------------------------------------------------------------- */
    int bRet = FALSE;
    if (y[0] != y[nPoints - 1] || y[0] != y[nMiddle] ||
        x[0] == x[nPoints - 1] || x[0] == x[nMiddle] ||
        (psATInfo->dfMaxErrorForward == 0.0 &&
         psATInfo->dfMaxErrorReverse == 0.0) ||
        nPoints <= 5)
    {
        bRet = psATInfo->pfnBaseTransformer(psATInfo->pBaseCBData, bDstToSrc,
                                            nPoints, x, y, z, panSuccess);
        goto end;
    }

    /* -------------------------------------------------------------------- */
    /*      Transform first, last and middle point.                         */
    /* -------------------------------------------------------------------- */
    x2[0] = x[0];
    y2[0] = y[0];
    z2[0] = z[0];
    x2[1] = x[nMiddle];
    y2[1] = y[nMiddle];
    z2[1] = z[nMiddle];
    x2[2] = x[nPoints - 1];
    y2[2] = y[nPoints - 1];
    z2[2] = z[nPoints - 1];

    bSuccess = psATInfo->pfnBaseTransformer(psATInfo->pBaseCBData, bDstToSrc, 3,
                                            x2, y2, z2, anSuccess2);
    if (!bSuccess || !anSuccess2[0] || !anSuccess2[1] || !anSuccess2[2])
    {
        bRet = psATInfo->pfnBaseTransformer(psATInfo->pBaseCBData, bDstToSrc,
                                            nPoints, x, y, z, panSuccess);
        goto end;
    }

    bRet = GDALApproxTransformInternal(pCBData, bDstToSrc, nPoints, x, y, z,
                                       panSuccess, x2, y2, z2);

end:
#ifdef DEBUG_APPROX_TRANSFORMER
    for (int i = 0; i < nPoints; i++)
        fprintf(stderr, "[%d] (%.10f,%.10f) %d\n", /*ok*/
                i, x[i], y[i], panSuccess[i]);
#endif

    return bRet;
}

/************************************************************************/
/*                  GDALDeserializeApproxTransformer()                  */
/************************************************************************/

static void *GDALDeserializeApproxTransformer(CPLXMLNode *psTree)

{
    double dfMaxErrorForward = 0.25;
    double dfMaxErrorReverse = 0.25;
    const char *pszMaxError = CPLGetXMLValue(psTree, "MaxError", nullptr);
    if (pszMaxError != nullptr)
    {
        dfMaxErrorForward = CPLAtof(pszMaxError);
        dfMaxErrorReverse = dfMaxErrorForward;
    }
    const char *pszMaxErrorForward =
        CPLGetXMLValue(psTree, "MaxErrorForward", nullptr);
    if (pszMaxErrorForward != nullptr)
    {
        dfMaxErrorForward = CPLAtof(pszMaxErrorForward);
    }
    const char *pszMaxErrorReverse =
        CPLGetXMLValue(psTree, "MaxErrorReverse", nullptr);
    if (pszMaxErrorReverse != nullptr)
    {
        dfMaxErrorReverse = CPLAtof(pszMaxErrorReverse);
    }

    GDALTransformerFunc pfnBaseTransform = nullptr;
    void *pBaseCBData = nullptr;

    CPLXMLNode *psContainer = CPLGetXMLNode(psTree, "BaseTransformer");

    if (psContainer != nullptr && psContainer->psChild != nullptr)
    {
        GDALDeserializeTransformer(psContainer->psChild, &pfnBaseTransform,
                                   &pBaseCBData);
    }

    if (pfnBaseTransform == nullptr)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Cannot get base transform for approx transformer.");
        return nullptr;
    }

    void *pApproxCBData = GDALCreateApproxTransformer2(
        pfnBaseTransform, pBaseCBData, dfMaxErrorForward, dfMaxErrorReverse);
    GDALApproxTransformerOwnsSubtransformer(pApproxCBData, TRUE);

    return pApproxCBData;
}

/************************************************************************/
/*                 GDALTransformLonLatToDestApproxTransformer()         */
/************************************************************************/

int GDALTransformLonLatToDestApproxTransformer(void *hTransformArg,
                                               double *pdfX, double *pdfY)
{
    GDALApproxTransformInfo *psInfo =
        static_cast<GDALApproxTransformInfo *>(hTransformArg);

    if (GDALIsTransformer(psInfo->pBaseCBData,
                          GDAL_GEN_IMG_TRANSFORMER_CLASS_NAME))
    {
        return GDALTransformLonLatToDestGenImgProjTransformer(
            psInfo->pBaseCBData, pdfX, pdfY);
    }
    return false;
}

/************************************************************************/
/*                       GDALApplyGeoTransform()                        */
/************************************************************************/

/**
 * Apply GeoTransform to x/y coordinate.
 *
 * Applies the following computation, converting a (pixel, line) coordinate
 * into a georeferenced (geo_x, geo_y) location.
 * \code{.c}
 *  *pdfGeoX = padfGeoTransform[0] + dfPixel * padfGeoTransform[1]
 *                                 + dfLine  * padfGeoTransform[2];
 *  *pdfGeoY = padfGeoTransform[3] + dfPixel * padfGeoTransform[4]
 *                                 + dfLine  * padfGeoTransform[5];
 * \endcode
 *
 * @param padfGeoTransform Six coefficient GeoTransform to apply.
 * @param dfPixel Input pixel position.
 * @param dfLine Input line position.
 * @param pdfGeoX output location where geo_x (easting/longitude)
 * location is placed.
 * @param pdfGeoY output location where geo_y (northing/latitude)
 * location is placed.
 */

void CPL_STDCALL GDALApplyGeoTransform(const double *padfGeoTransform,
                                       double dfPixel, double dfLine,
                                       double *pdfGeoX, double *pdfGeoY)
{
    *pdfGeoX = padfGeoTransform[0] + dfPixel * padfGeoTransform[1] +
               dfLine * padfGeoTransform[2];
    *pdfGeoY = padfGeoTransform[3] + dfPixel * padfGeoTransform[4] +
               dfLine * padfGeoTransform[5];
}

/************************************************************************/
/*                        GDALInvGeoTransform()                         */
/************************************************************************/

/**
 * Invert Geotransform.
 *
 * This function will invert a standard 3x2 set of GeoTransform coefficients.
 * This converts the equation from being pixel to geo to being geo to pixel.
 *
 * @param gt_in Input geotransform (six doubles - unaltered).
 * @param gt_out Output geotransform (six doubles - updated).
 *
 * @return TRUE on success or FALSE if the equation is uninvertable.
 */

int CPL_STDCALL GDALInvGeoTransform(const double *gt_in, double *gt_out)

{
    // Special case - no rotation - to avoid computing determinate
    // and potential precision issues.
    if (gt_in[2] == 0.0 && gt_in[4] == 0.0 && gt_in[1] != 0.0 &&
        gt_in[5] != 0.0)
    {
        /*X = gt_in[0] + x * gt_in[1]
          Y = gt_in[3] + y * gt_in[5]
          -->
          x = -gt_in[0] / gt_in[1] + (1 / gt_in[1]) * X
          y = -gt_in[3] / gt_in[5] + (1 / gt_in[5]) * Y
        */
        gt_out[0] = -gt_in[0] / gt_in[1];
        gt_out[1] = 1.0 / gt_in[1];
        gt_out[2] = 0.0;
        gt_out[3] = -gt_in[3] / gt_in[5];
        gt_out[4] = 0.0;
        gt_out[5] = 1.0 / gt_in[5];
        return 1;
    }

    // Assume a 3rd row that is [1 0 0].

    // Compute determinate.

    const double det = gt_in[1] * gt_in[5] - gt_in[2] * gt_in[4];
    const double magnitude = std::max(std::max(fabs(gt_in[1]), fabs(gt_in[2])),
                                      std::max(fabs(gt_in[4]), fabs(gt_in[5])));

    if (fabs(det) <= 1e-10 * magnitude * magnitude)
        return 0;

    const double inv_det = 1.0 / det;

    // Compute adjoint, and divide by determinate.

    gt_out[1] = gt_in[5] * inv_det;
    gt_out[4] = -gt_in[4] * inv_det;

    gt_out[2] = -gt_in[2] * inv_det;
    gt_out[5] = gt_in[1] * inv_det;

    gt_out[0] = (gt_in[2] * gt_in[3] - gt_in[0] * gt_in[5]) * inv_det;
    gt_out[3] = (-gt_in[1] * gt_in[3] + gt_in[0] * gt_in[4]) * inv_det;

    return 1;
}

/************************************************************************/
/*                      GDALSerializeTransformer()                      */
/************************************************************************/

CPLXMLNode *GDALSerializeTransformer(GDALTransformerFunc /* pfnFunc */,
                                     void *pTransformArg)
{
    VALIDATE_POINTER1(pTransformArg, "GDALSerializeTransformer", nullptr);

    GDALTransformerInfo *psInfo =
        static_cast<GDALTransformerInfo *>(pTransformArg);

    if (psInfo == nullptr || memcmp(psInfo->abySignature, GDAL_GTI2_SIGNATURE,
                                    strlen(GDAL_GTI2_SIGNATURE)) != 0)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Attempt to serialize non-GTI2 transformer.");
        return nullptr;
    }
    else if (psInfo->pfnSerialize == nullptr)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "No serialization function available for this transformer.");
        return nullptr;
    }

    return psInfo->pfnSerialize(pTransformArg);
}

/************************************************************************/
/*                  GDALRegisterTransformDeserializer()                 */
/************************************************************************/

static CPLList *psListDeserializer = nullptr;
static CPLMutex *hDeserializerMutex = nullptr;

typedef struct
{
    char *pszTransformName;
    GDALTransformerFunc pfnTransformerFunc;
    GDALTransformDeserializeFunc pfnDeserializeFunc;
} TransformDeserializerInfo;

void *GDALRegisterTransformDeserializer(
    const char *pszTransformName, GDALTransformerFunc pfnTransformerFunc,
    GDALTransformDeserializeFunc pfnDeserializeFunc)
{
    TransformDeserializerInfo *psInfo =
        static_cast<TransformDeserializerInfo *>(
            CPLMalloc(sizeof(TransformDeserializerInfo)));
    psInfo->pszTransformName = CPLStrdup(pszTransformName);
    psInfo->pfnTransformerFunc = pfnTransformerFunc;
    psInfo->pfnDeserializeFunc = pfnDeserializeFunc;

    CPLMutexHolderD(&hDeserializerMutex);
    psListDeserializer = CPLListInsert(psListDeserializer, psInfo, 0);

    return psInfo;
}

/************************************************************************/
/*                GDALUnregisterTransformDeserializer()                 */
/************************************************************************/

void GDALUnregisterTransformDeserializer(void *pData)
{
    CPLMutexHolderD(&hDeserializerMutex);
    CPLList *psList = psListDeserializer;
    CPLList *psLast = nullptr;
    while (psList)
    {
        if (psList->pData == pData)
        {
            TransformDeserializerInfo *psInfo =
                static_cast<TransformDeserializerInfo *>(pData);
            CPLFree(psInfo->pszTransformName);
            CPLFree(pData);
            if (psLast)
                psLast->psNext = psList->psNext;
            else
                psListDeserializer = nullptr;
            CPLFree(psList);
            break;
        }
        psLast = psList;
        psList = psList->psNext;
    }
}

/************************************************************************/
/*                GDALUnregisterTransformDeserializer()                 */
/************************************************************************/

void GDALCleanupTransformDeserializerMutex()
{
    if (hDeserializerMutex != nullptr)
    {
        CPLDestroyMutex(hDeserializerMutex);
        hDeserializerMutex = nullptr;
    }
}

/************************************************************************/
/*                     GDALDeserializeTransformer()                     */
/************************************************************************/

CPLErr GDALDeserializeTransformer(CPLXMLNode *psTree,
                                  GDALTransformerFunc *ppfnFunc,
                                  void **ppTransformArg)

{
    *ppfnFunc = nullptr;
    *ppTransformArg = nullptr;

    CPLErrorReset();

    if (psTree == nullptr || psTree->eType != CXT_Element)
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Malformed element in GDALDeserializeTransformer");
    else if (EQUAL(psTree->pszValue, "GenImgProjTransformer"))
    {
        *ppfnFunc = GDALGenImgProjTransform;
        *ppTransformArg = GDALDeserializeGenImgProjTransformer(psTree);
    }
    else if (EQUAL(psTree->pszValue, "ReprojectionTransformer"))
    {
        *ppfnFunc = GDALReprojectionTransform;
        *ppTransformArg = GDALDeserializeReprojectionTransformer(psTree);
    }
    else if (EQUAL(psTree->pszValue, "GCPTransformer"))
    {
        *ppfnFunc = GDALGCPTransform;
        *ppTransformArg = GDALDeserializeGCPTransformer(psTree);
    }
    else if (EQUAL(psTree->pszValue, "TPSTransformer"))
    {
        *ppfnFunc = GDALTPSTransform;
        *ppTransformArg = GDALDeserializeTPSTransformer(psTree);
    }
    else if (EQUAL(psTree->pszValue, "GeoLocTransformer"))
    {
        *ppfnFunc = GDALGeoLocTransform;
        *ppTransformArg = GDALDeserializeGeoLocTransformer(psTree);
    }
    else if (EQUAL(psTree->pszValue, "RPCTransformer"))
    {
        *ppfnFunc = GDALRPCTransform;
        *ppTransformArg = GDALDeserializeRPCTransformer(psTree);
    }
    else if (EQUAL(psTree->pszValue, "ApproxTransformer"))
    {
        *ppfnFunc = GDALApproxTransform;
        *ppTransformArg = GDALDeserializeApproxTransformer(psTree);
    }
    else if (EQUAL(psTree->pszValue, "HomographyTransformer"))
    {
        *ppfnFunc = GDALHomographyTransform;
        *ppTransformArg = GDALDeserializeHomographyTransformer(psTree);
    }
    else
    {
        GDALTransformDeserializeFunc pfnDeserializeFunc = nullptr;
        {
            CPLMutexHolderD(&hDeserializerMutex);
            CPLList *psList = psListDeserializer;
            while (psList)
            {
                TransformDeserializerInfo *psInfo =
                    static_cast<TransformDeserializerInfo *>(psList->pData);
                if (strcmp(psInfo->pszTransformName, psTree->pszValue) == 0)
                {
                    *ppfnFunc = psInfo->pfnTransformerFunc;
                    pfnDeserializeFunc = psInfo->pfnDeserializeFunc;
                    break;
                }
                psList = psList->psNext;
            }
        }

        if (pfnDeserializeFunc != nullptr)
        {
            *ppTransformArg = pfnDeserializeFunc(psTree);
        }
        else
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Unrecognized element '%s' GDALDeserializeTransformer",
                     psTree->pszValue);
        }
    }

    return CPLGetLastErrorType();
}

/************************************************************************/
/*                       GDALDestroyTransformer()                       */
/************************************************************************/

void GDALDestroyTransformer(void *pTransformArg)

{
    if (pTransformArg == nullptr)
        return;

    GDALTransformerInfo *psInfo =
        static_cast<GDALTransformerInfo *>(pTransformArg);

    if (memcmp(psInfo->abySignature, GDAL_GTI2_SIGNATURE,
               strlen(GDAL_GTI2_SIGNATURE)) != 0)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Attempt to destroy non-GTI2 transformer.");
        return;
    }

    psInfo->pfnCleanup(pTransformArg);
}

/************************************************************************/
/*                         GDALUseTransformer()                         */
/************************************************************************/

int GDALUseTransformer(void *pTransformArg, int bDstToSrc, int nPointCount,
                       double *x, double *y, double *z, int *panSuccess)
{
    GDALTransformerInfo *psInfo =
        static_cast<GDALTransformerInfo *>(pTransformArg);

    if (psInfo == nullptr || memcmp(psInfo->abySignature, GDAL_GTI2_SIGNATURE,
                                    strlen(GDAL_GTI2_SIGNATURE)) != 0)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Attempt to use non-GTI2 transformer.");
        return FALSE;
    }

    return psInfo->pfnTransform(pTransformArg, bDstToSrc, nPointCount, x, y, z,
                                panSuccess);
}

/************************************************************************/
/*                        GDALCloneTransformer()                        */
/************************************************************************/

void *GDALCloneTransformer(void *pTransformArg)
{
    GDALTransformerInfo *psInfo =
        static_cast<GDALTransformerInfo *>(pTransformArg);

    if (psInfo == nullptr || memcmp(psInfo->abySignature, GDAL_GTI2_SIGNATURE,
                                    strlen(GDAL_GTI2_SIGNATURE)) != 0)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Attempt to clone non-GTI2 transformer.");
        return nullptr;
    }

    if (psInfo->pfnCreateSimilar != nullptr)
    {
        return psInfo->pfnCreateSimilar(psInfo, 1.0, 1.0);
    }

    if (psInfo->pfnSerialize == nullptr)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "No serialization function available for this transformer.");
        return nullptr;
    }

    CPLXMLNode *pSerialized = psInfo->pfnSerialize(pTransformArg);
    if (pSerialized == nullptr)
        return nullptr;
    GDALTransformerFunc pfnTransformer = nullptr;
    void *pClonedTransformArg = nullptr;
    if (GDALDeserializeTransformer(pSerialized, &pfnTransformer,
                                   &pClonedTransformArg) != CE_None)
    {
        CPLDestroyXMLNode(pSerialized);
        CPLFree(pClonedTransformArg);
        return nullptr;
    }

    CPLDestroyXMLNode(pSerialized);
    return pClonedTransformArg;
}

/************************************************************************/
/*                   GDALCreateSimilarTransformer()                     */
/************************************************************************/

void *GDALCreateSimilarTransformer(void *pTransformArg, double dfRatioX,
                                   double dfRatioY)
{
    GDALTransformerInfo *psInfo =
        static_cast<GDALTransformerInfo *>(pTransformArg);

    if (psInfo == nullptr || memcmp(psInfo->abySignature, GDAL_GTI2_SIGNATURE,
                                    strlen(GDAL_GTI2_SIGNATURE)) != 0)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Attempt to call CreateSimilar on a non-GTI2 transformer.");
        return nullptr;
    }

    if (psInfo->pfnCreateSimilar == nullptr)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "No CreateSimilar function available for this transformer.");
        return nullptr;
    }

    return psInfo->pfnCreateSimilar(psInfo, dfRatioX, dfRatioY);
}

/************************************************************************/
/*                      GetGenImgProjTransformInfo()                    */
/************************************************************************/

static GDALTransformerInfo *GetGenImgProjTransformInfo(const char *pszFunc,
                                                       void *pTransformArg)
{
    GDALTransformerInfo *psInfo =
        static_cast<GDALTransformerInfo *>(pTransformArg);

    if (psInfo == nullptr || memcmp(psInfo->abySignature, GDAL_GTI2_SIGNATURE,
                                    strlen(GDAL_GTI2_SIGNATURE)) != 0)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Attempt to call %s on "
                 "a non-GTI2 transformer.",
                 pszFunc);
        return nullptr;
    }

    if (EQUAL(psInfo->pszClassName, GDAL_APPROX_TRANSFORMER_CLASS_NAME))
    {
        GDALApproxTransformInfo *psATInfo =
            static_cast<GDALApproxTransformInfo *>(pTransformArg);
        psInfo = static_cast<GDALTransformerInfo *>(psATInfo->pBaseCBData);

        if (psInfo == nullptr ||
            memcmp(psInfo->abySignature, GDAL_GTI2_SIGNATURE,
                   strlen(GDAL_GTI2_SIGNATURE)) != 0)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Attempt to call %s on "
                     "a non-GTI2 transformer.",
                     pszFunc);
            return nullptr;
        }
    }

    if (EQUAL(psInfo->pszClassName, GDAL_GEN_IMG_TRANSFORMER_CLASS_NAME))
    {
        return psInfo;
    }

    return nullptr;
}

/************************************************************************/
/*                 GDALSetTransformerDstGeoTransform()                  */
/************************************************************************/

/**
 * Set ApproxTransformer or GenImgProj output geotransform.
 *
 * This is a layer above GDALSetGenImgProjTransformerDstGeoTransform() that
 * checks that the passed hTransformArg is compatible.
 *
 * Normally the "destination geotransform", or transformation between
 * georeferenced output coordinates and pixel/line coordinates on the
 * destination file is extracted from the destination file by
 * GDALCreateGenImgProjTransformer() and stored in the GenImgProj private
 * info.  However, sometimes it is inconvenient to have an output file
 * handle with appropriate geotransform information when creating the
 * transformation.  For these cases, this function can be used to apply
 * the destination geotransform.
 *
 * @param pTransformArg the handle to update.
 * @param padfGeoTransform the destination geotransform to apply (six doubles).
 */

void GDALSetTransformerDstGeoTransform(void *pTransformArg,
                                       const double *padfGeoTransform)
{
    VALIDATE_POINTER0(pTransformArg, "GDALSetTransformerDstGeoTransform");

    GDALTransformerInfo *psInfo = GetGenImgProjTransformInfo(
        "GDALSetTransformerDstGeoTransform", pTransformArg);
    if (psInfo)
    {
        GDALSetGenImgProjTransformerDstGeoTransform(psInfo, padfGeoTransform);
    }
}

/************************************************************************/
/*                 GDALGetTransformerDstGeoTransform()                  */
/************************************************************************/

/**
 * Get ApproxTransformer or GenImgProj output geotransform.
 *
 * @param pTransformArg transformer handle.
 * @param padfGeoTransform (output) the destination geotransform to return (six
 * doubles).
 */

void GDALGetTransformerDstGeoTransform(void *pTransformArg,
                                       double *padfGeoTransform)
{
    VALIDATE_POINTER0(pTransformArg, "GDALGetTransformerDstGeoTransform");

    GDALTransformerInfo *psInfo = GetGenImgProjTransformInfo(
        "GDALGetTransformerDstGeoTransform", pTransformArg);
    if (psInfo)
    {
        GDALGenImgProjTransformInfo *psGenImgProjInfo =
            reinterpret_cast<GDALGenImgProjTransformInfo *>(psInfo);

        memcpy(padfGeoTransform, psGenImgProjInfo->sDstParams.adfGeoTransform,
               sizeof(double) * 6);
    }
}

/************************************************************************/
/*            GDALTransformIsTranslationOnPixelBoundaries()             */
/************************************************************************/

bool GDALTransformIsTranslationOnPixelBoundaries(GDALTransformerFunc,
                                                 void *pTransformerArg)
{
    if (GDALIsTransformer(pTransformerArg, GDAL_APPROX_TRANSFORMER_CLASS_NAME))
    {
        const auto *pApproxInfo =
            static_cast<const GDALApproxTransformInfo *>(pTransformerArg);
        pTransformerArg = pApproxInfo->pBaseCBData;
    }
    if (GDALIsTransformer(pTransformerArg, GDAL_GEN_IMG_TRANSFORMER_CLASS_NAME))
    {
        const auto *pGenImgpProjInfo =
            static_cast<GDALGenImgProjTransformInfo *>(pTransformerArg);
        const auto IsCloseToInteger = [](double dfVal)
        { return std::fabs(dfVal - std::round(dfVal)) <= 1e-6; };
        return pGenImgpProjInfo->sSrcParams.pTransformArg == nullptr &&
               pGenImgpProjInfo->sDstParams.pTransformArg == nullptr &&
               pGenImgpProjInfo->pReproject == nullptr &&
               pGenImgpProjInfo->sSrcParams.adfGeoTransform[1] ==
                   pGenImgpProjInfo->sDstParams.adfGeoTransform[1] &&
               pGenImgpProjInfo->sSrcParams.adfGeoTransform[5] ==
                   pGenImgpProjInfo->sDstParams.adfGeoTransform[5] &&
               pGenImgpProjInfo->sSrcParams.adfGeoTransform[2] ==
                   pGenImgpProjInfo->sDstParams.adfGeoTransform[2] &&
               pGenImgpProjInfo->sSrcParams.adfGeoTransform[4] ==
                   pGenImgpProjInfo->sDstParams.adfGeoTransform[4] &&
               // Check that the georeferenced origin of the destination
               // geotransform is close to be an integer value when transformed
               // to source image coordinates
               IsCloseToInteger(
                   pGenImgpProjInfo->sSrcParams.adfInvGeoTransform[0] +
                   pGenImgpProjInfo->sDstParams.adfGeoTransform[0] *
                       pGenImgpProjInfo->sSrcParams.adfInvGeoTransform[1] +
                   pGenImgpProjInfo->sDstParams.adfGeoTransform[3] *
                       pGenImgpProjInfo->sSrcParams.adfInvGeoTransform[2]) &&
               IsCloseToInteger(
                   pGenImgpProjInfo->sSrcParams.adfInvGeoTransform[3] +
                   pGenImgpProjInfo->sDstParams.adfGeoTransform[0] *
                       pGenImgpProjInfo->sSrcParams.adfInvGeoTransform[4] +
                   pGenImgpProjInfo->sDstParams.adfGeoTransform[3] *
                       pGenImgpProjInfo->sSrcParams.adfInvGeoTransform[5]);
    }
    return false;
}

/************************************************************************/
/*                   GDALTransformIsAffineNoRotation()                  */
/************************************************************************/

bool GDALTransformIsAffineNoRotation(GDALTransformerFunc, void *pTransformerArg)
{
    if (GDALIsTransformer(pTransformerArg, GDAL_APPROX_TRANSFORMER_CLASS_NAME))
    {
        const auto *pApproxInfo =
            static_cast<const GDALApproxTransformInfo *>(pTransformerArg);
        pTransformerArg = pApproxInfo->pBaseCBData;
    }
    if (GDALIsTransformer(pTransformerArg, GDAL_GEN_IMG_TRANSFORMER_CLASS_NAME))
    {
        const auto *pGenImgpProjInfo =
            static_cast<GDALGenImgProjTransformInfo *>(pTransformerArg);
        return pGenImgpProjInfo->sSrcParams.pTransformArg == nullptr &&
               pGenImgpProjInfo->sDstParams.pTransformArg == nullptr &&
               pGenImgpProjInfo->pReproject == nullptr &&
               pGenImgpProjInfo->sSrcParams.adfGeoTransform[2] == 0 &&
               pGenImgpProjInfo->sSrcParams.adfGeoTransform[4] == 0 &&
               pGenImgpProjInfo->sDstParams.adfGeoTransform[2] == 0 &&
               pGenImgpProjInfo->sDstParams.adfGeoTransform[4] == 0;
    }
    return false;
}

/************************************************************************/
/*                        GDALTransformHasFastClone()                   */
/************************************************************************/

/** Returns whether GDALCloneTransformer() on this transformer is
 * "fast"
 * Counter-examples are GCPs or TPSs transformers.
 */
bool GDALTransformHasFastClone(void *pTransformerArg)
{
    if (GDALIsTransformer(pTransformerArg, GDAL_APPROX_TRANSFORMER_CLASS_NAME))
    {
        const auto *pApproxInfo =
            static_cast<const GDALApproxTransformInfo *>(pTransformerArg);
        pTransformerArg = pApproxInfo->pBaseCBData;
        // Fallback to next lines
    }

    if (GDALIsTransformer(pTransformerArg, GDAL_GEN_IMG_TRANSFORMER_CLASS_NAME))
    {
        const auto *pGenImgpProjInfo =
            static_cast<GDALGenImgProjTransformInfo *>(pTransformerArg);
        return (pGenImgpProjInfo->sSrcParams.pTransformArg == nullptr ||
                GDALTransformHasFastClone(
                    pGenImgpProjInfo->sSrcParams.pTransformArg)) &&
               (pGenImgpProjInfo->sDstParams.pTransformArg == nullptr ||
                GDALTransformHasFastClone(
                    pGenImgpProjInfo->sDstParams.pTransformArg));
    }
    else if (GDALIsTransformer(pTransformerArg,
                               GDAL_RPC_TRANSFORMER_CLASS_NAME))
    {
        return true;
    }
    else
    {
        return false;
    }
}
