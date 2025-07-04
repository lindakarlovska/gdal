/******************************************************************************
 *
 * Project:  MVT Translator
 * Purpose:  Coordinate and envelope transformations related to MVT.
 * Author:   Even Rouault, Even Rouault <even dot rouault at spatialys dot com>
 *           Linda Karlovska <linda dot karlovska at seznam dot cz>
 *
 ******************************************************************************
 * Copyright (c) 2018, Even Rouault <even dot rouault at spatialys dot com>
 * Copyright (c) 2025, Linda Karlovska <linda dot karlovska at seznam dot cz>
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#include "geoutils.h"

#include "cpl_conv.h"
#include "ogr_geos.h"
#include "ogr_p.h"

#include <cmath>
#include <algorithm>
#include <iostream>

/************************************************************************/
/*                    InitWebMercatorTilingScheme()                     */
/************************************************************************/
void MVTGeoUtils::InitWebMercatorTilingScheme(OGRSpatialReference *poSRS,
                                              double &dfTopX, double &dfTopY,
                                              double &dfTileDim0)
{
    std::cout << "InitWebMercatorTilingScheme" << std::endl;

    poSRS->SetFromUserInput(kSRSWebMercator);
    dfTopX = -kmMAX_GM;
    dfTopY = kmMAX_GM;
    dfTileDim0 = 2 * kmMAX_GM;
}

/************************************************************************/
/*                    InitUserDefinedTilingScheme()                     */
/************************************************************************/
bool MVTGeoUtils::InitUserDefinedTilingScheme(const char *pszTilingScheme,
                                              OGRSpatialReference *poSRS,
                                              double &dfTopX, double &dfTopY,
                                              double &dfTileDim0,
                                              int &nTileMatrixWidth0,
                                              int &nTileMatrixHeight0)
{
    std::cout << "InitUserDefinedTilingScheme" << std::endl;

    CPLStringList aoList(CSLTokenizeString2(pszTilingScheme, ",", 0));
    if (aoList.Count() >= 4)
    {
        poSRS->SetFromUserInput(aoList[0]);
        dfTopX = CPLAtof(aoList[1]);
        dfTopY = CPLAtof(aoList[2]);
        dfTileDim0 = CPLAtof(aoList[3]);

        if (aoList.Count() == 6)
        {
            nTileMatrixWidth0 = std::max(1, atoi(aoList[4]));
            nTileMatrixHeight0 = std::max(1, atoi(aoList[5]));
        }
        else if (dfTopX == -180 && dfTileDim0 == 180)
        {
            // Assumes WorldCRS84Quad
            nTileMatrixWidth0 = 2;
            nTileMatrixHeight0 = 1;
        }
        else
        {
            // default fallback
            nTileMatrixWidth0 = 1;
            nTileMatrixHeight0 = 1;
        }
        return true;
    }
    else
    {
        CPLError(
            CE_Failure, CPLE_AppDefined,
            "Wrong format for TILING_SCHEME. "
            "Expecting "
            "EPSG:XXXX,tile_origin_x,tile_origin_y,tile_dim[,width,height]");
        return false;
    }
}

/************************************************************************/
/*                     SphericalMercatorToLongLat()                     */
/************************************************************************/
void MVTGeoUtils::SphericalMercatorToLongLat(double *pdfX, double *pdfY)
{
    std::cout << "SphericalMercatorToLongLat" << std::endl;
    double lng = *pdfX / kmSPHERICAL_RADIUS / M_PI * 180;
    double lat =
        2 * (atan(exp(*pdfY / kmSPHERICAL_RADIUS)) - M_PI / 4) / M_PI * 180;
    *pdfX = lng;
    *pdfY = lat;
}

/************************************************************************/
/*                     LongLatToSphericalMercator()                     */
/************************************************************************/
void MVTGeoUtils::LongLatToSphericalMercator(double *pdfX, double *pdfY)
{
    std::cout << "LongLatToSphericalMercator" << std::endl;

    double X = kmSPHERICAL_RADIUS * (*pdfX) / 180 * M_PI;
    double Y =
        kmSPHERICAL_RADIUS * log(tan(M_PI / 4 + 0.5 * (*pdfY) / 180 * M_PI));
    *pdfX = X;
    *pdfY = Y;
}

/************************************************************************/
/*                       ConvertFromWGS84()                             */
/************************************************************************/
void MVTGeoUtils::ConvertFromWGS84(OGRSpatialReference *poTargetSRS,
                                   double &dfX0, double &dfY0, double &dfX1,
                                   double &dfY1)
{
    std::cout << "ConvertFromWGS84" << std::endl;
    OGRSpatialReference oSRS_EPSG3857;
    oSRS_EPSG3857.SetFromUserInput(kSRSWebMercator);

    if (poTargetSRS->IsSame(&oSRS_EPSG3857))
    {
        LongLatToSphericalMercator(&dfX0, &dfY0);
        LongLatToSphericalMercator(&dfX1, &dfY1);
    }
    else
    {
        OGRSpatialReference oSRS_EPSG4326;
        oSRS_EPSG4326.SetFromUserInput(kSRSWGS84);
        oSRS_EPSG4326.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        OGRCoordinateTransformation *poCT =
            OGRCreateCoordinateTransformation(&oSRS_EPSG4326, poTargetSRS);
        if (poCT)
        {
            poCT->Transform(1, &dfX0, &dfY0);
            poCT->Transform(1, &dfX1, &dfY1);
            delete poCT;
        }
    }
}

/************************************************************************/
/*                       ComputeStandardEnvelope()                      */
/************************************************************************/
void MVTGeoUtils::ComputeStandardEnvelope(const OGREnvelope &oEnvelope,
                                          OGREnvelope &oTransformedEnvelope)
{
    std::cout << "ComputeStandardEnvelope" << std::endl;

    oTransformedEnvelope = oEnvelope;

    SphericalMercatorToLongLat(&(oTransformedEnvelope.MinX),
                               &(oTransformedEnvelope.MinY));
    SphericalMercatorToLongLat(&(oTransformedEnvelope.MaxX),
                               &(oTransformedEnvelope.MaxY));

    oTransformedEnvelope.MinY = std::max(-85.0, oTransformedEnvelope.MinY);
    oTransformedEnvelope.MaxY = std::min(85.0, oTransformedEnvelope.MaxY);
}

/************************************************************************/
/*                  TransformEnvelopeToWGS84()                          */
/************************************************************************/
void MVTGeoUtils::TransformEnvelopeToWGS84(const OGREnvelope &oEnvelope,
                                           OGRSpatialReference *poSRS,
                                           OGREnvelope &oTransformedEnvelope)
{
    std::cout << "TransformEnvelopeToWGS84" << std::endl;

    oTransformedEnvelope = oEnvelope;

    OGRSpatialReference oSRS_EPSG4326;
    oSRS_EPSG4326.SetFromUserInput(kSRSWGS84);
    oSRS_EPSG4326.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation *poCT =
        OGRCreateCoordinateTransformation(poSRS, &oSRS_EPSG4326);

    if (poCT)
    {
        OGRPoint oPoint1(oTransformedEnvelope.MinX, oTransformedEnvelope.MinY);
        OGRPoint oPoint2(oTransformedEnvelope.MinX, oTransformedEnvelope.MaxY);
        OGRPoint oPoint3(oTransformedEnvelope.MaxX, oTransformedEnvelope.MaxY);
        OGRPoint oPoint4(oTransformedEnvelope.MaxX, oTransformedEnvelope.MinY);

        oPoint1.transform(poCT);
        oPoint2.transform(poCT);
        oPoint3.transform(poCT);
        oPoint4.transform(poCT);

        oTransformedEnvelope.MinX = std::min(
            {oPoint1.getX(), oPoint2.getX(), oPoint3.getX(), oPoint4.getX()});
        oTransformedEnvelope.MinY = std::min(
            {oPoint1.getY(), oPoint2.getY(), oPoint3.getY(), oPoint4.getY()});
        oTransformedEnvelope.MaxX = std::max(
            {oPoint1.getX(), oPoint2.getX(), oPoint3.getX(), oPoint4.getX()});
        oTransformedEnvelope.MaxY = std::max(
            {oPoint1.getY(), oPoint2.getY(), oPoint3.getY(), oPoint4.getY()});

        delete poCT;
    }
}

/************************************************************************/
/*                           ComputeCenter()                            */
/************************************************************************/

CPLString MVTGeoUtils::ComputeCenter(const OGREnvelope &oEnvelope, int nMinZoom)
{
    std::cout << "ComputeCenter" << std::endl;

    const double dfCenterX = (oEnvelope.MinX + oEnvelope.MaxX) / 2.0;
    const double dfCenterY = (oEnvelope.MinY + oEnvelope.MaxY) / 2.0;

    return CPLSPrintf("%.7f,%.7f,%d", dfCenterX, dfCenterY, nMinZoom);
}

/************************************************************************/
/*                          ComputeBounds()                    		    */
/************************************************************************/
CPLString MVTGeoUtils::ComputeBounds(const OGREnvelope &oEnvelope)
{
    std::cout << "ComputeBounds" << std::endl;

    return CPLSPrintf("%.7f,%.7f,%.7f,%.7f", oEnvelope.MinX, oEnvelope.MinY,
                      oEnvelope.MaxX, oEnvelope.MaxY);
}

/************************************************************************/
/*                          ValidateMinMaxZoom()                        */
/************************************************************************/
bool MVTGeoUtils::ValidateMinMaxZoom(int nMinZoom, int nMaxZoom)
{
    if (nMinZoom < 0 || nMinZoom > 22)
    {
        CPLError(CE_Failure, CPLE_AppDefined, "Invalid MINZOOM");
        return false;
    }
    if (nMaxZoom < 0 || nMaxZoom > 22)
    {
        CPLError(CE_Failure, CPLE_AppDefined, "Invalid MAXZOOM");
        return false;
    }
    if (nMaxZoom < nMinZoom)
    {
        CPLError(CE_Failure, CPLE_AppDefined, "Invalid MAXZOOM < MINZOOM");
        return false;
    }
    return true;
}
