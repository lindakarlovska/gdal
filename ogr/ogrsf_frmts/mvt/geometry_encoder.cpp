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

#include "geometry_encoder.h"
#include "mvtutils.h"

#include <iostream>

/************************************************************************/
/*               GeometryEncoder::GeometryEncoder()                     */
/************************************************************************/
GeometryEncoder::GeometryEncoder(double dfTileDim, double dfTopX, double dfTopY,
                                 int nExtent)
    : m_dfTileDim(dfTileDim), m_dfTopX(dfTopX), m_dfTopY(dfTopY),
      m_nExtent(nExtent)
{
}

/************************************************************************/
/*                          	Encode()                                */
/************************************************************************/
bool GeometryEncoder::Encode(const OGRGeometry *poOriginalGeom,
                             const OGRGeometry *poGeomToEncode,
                             MVTTileLayerFeature *poGPBFeature,
                             double &dfAreaOrLength) const
{
    std::cout << "EncodeGeometry" << std::endl;

    bool bGeomOK = false;

    // Get the original geometry type for decisions
    auto eGeomType = wkbFlatten(poOriginalGeom->getGeometryType());
    // Get the type of the geometry to encode
    auto eGeomToEncodeType = wkbFlatten(poGeomToEncode->getGeometryType());

    if (eGeomType == wkbPoint || eGeomType == wkbMultiPoint)
    {
        if (eGeomToEncodeType == wkbPoint)
        {
            const OGRPoint *poPoint = poGeomToEncode->toPoint();
            bGeomOK = HandlePointGeometry(poPoint, poGPBFeature);
        }
        else if (eGeomToEncodeType == wkbMultiPoint ||
                 eGeomToEncodeType == wkbGeometryCollection)
        {
            const OGRGeometryCollection *poGC =
                poGeomToEncode->toGeometryCollection();
            bGeomOK = HandleMultiPointOrCollection(poGC, poGPBFeature);
        }
    }
    else if (eGeomType == wkbLineString || eGeomType == wkbMultiLineString)
    {
        if (eGeomToEncodeType == wkbLineString)
        {
            const OGRLineString *poLS = poGeomToEncode->toLineString();
            bGeomOK =
                HandleLineStringGeometry(poLS, poGPBFeature, dfAreaOrLength);
        }
        else if (eGeomToEncodeType == wkbMultiLineString ||
                 eGeomToEncodeType == wkbGeometryCollection)
        {
            const OGRGeometryCollection *poGC =
                poGeomToEncode->toGeometryCollection();
            bGeomOK = HandleMultiLineStringOrCollection(poGC, poGPBFeature,
                                                        dfAreaOrLength);
        }
    }
    else if (eGeomType == wkbPolygon || eGeomType == wkbMultiPolygon)
    {
        if (eGeomToEncodeType == wkbPolygon)
        {
            const OGRPolygon *poPoly = poGeomToEncode->toPolygon();
            bGeomOK =
                HandlePolygonGeometry(poPoly, poGPBFeature, dfAreaOrLength);
        }
        else if (eGeomToEncodeType == wkbMultiPolygon ||
                 eGeomToEncodeType == wkbGeometryCollection)
        {
            const OGRGeometryCollection *poGC =
                poGeomToEncode->toGeometryCollection();
            bGeomOK = HandleMultiPolygonOrCollection(poGC, poGPBFeature,
                                                     dfAreaOrLength);
        }
    }

    return bGeomOK;
}

/************************************************************************/
/*                        ConvertToTileCoords()                     	*/
/************************************************************************/
void GeometryEncoder::ConvertToTileCoords(double dfX, double dfY, int &nX,
                                          int &nY, double dfTopX, double dfTopY,
                                          double dfTileDim) const
{
    if (dfTileDim == 0)
    {
        nX = static_cast<int>(dfX);
        nY = static_cast<int>(dfY);
    }
    else
    {
        nX = static_cast<int>(
            std::round((dfX - dfTopX) * m_nExtent / dfTileDim));
        nY = static_cast<int>(
            std::round((dfTopY - dfY) * m_nExtent / dfTileDim));
    }
}

/************************************************************************/
/*                        ConvertToTileCoords()                     	*/
/************************************************************************/
void GeometryEncoder::ConvertToTileCoords(double dfX, double dfY, int &nX,
                                          int &nY) const
{
    ConvertToTileCoords(dfX, dfY, nX, nY, m_dfTopX, m_dfTopY, m_dfTileDim);
}

/************************************************************************/
/*                          HandlePointGeometry()                       */
/************************************************************************/
bool GeometryEncoder::HandlePointGeometry(
    const OGRPoint *poPoint, MVTTileLayerFeature *poGPBFeature) const
{
    std::cout << "HandlePointGeometry" << std::endl;

    int nX, nY;
    double dfX = poPoint->getX();
    double dfY = poPoint->getY();

    // Convert to tile coordinates
    ConvertToTileCoords(dfX, dfY, nX, nY);

    // Add MOVETO command and encoded point coordinates
    poGPBFeature->addGeometry(GetCmdCountCombined(knCMD_MOVETO, 1));
    poGPBFeature->addGeometry(EncodeSInt(nX));
    poGPBFeature->addGeometry(EncodeSInt(nY));

    return true;
}

/************************************************************************/
/*                          HandlePointCollection()                     */
/************************************************************************/
bool GeometryEncoder::HandleMultiPointOrCollection(
    const OGRGeometryCollection *poGC, MVTTileLayerFeature *poGPBFeature) const
{
    std::cout << "HandlePointCollection" << std::endl;

    std::set<std::pair<int, int>> oSetUniqueCoords;
    poGPBFeature->addGeometry(
        GetCmdCountCombined(knCMD_MOVETO, 0));  // To be modified later
    int nLastX = 0;
    int nLastY = 0;

    for (auto &&poSubGeom : poGC)
    {
        if (wkbFlatten(poSubGeom->getGeometryType()) == wkbPoint)
        {
            const OGRPoint *poPoint = poSubGeom->toPoint();
            int nX, nY;
            double dfX = poPoint->getX();
            double dfY = poPoint->getY();

            // Convert to tile coordinates
            ConvertToTileCoords(dfX, dfY, nX, nY);

            // Check if the coordinate is unique
            if (oSetUniqueCoords.find(std::pair<int, int>(nX, nY)) ==
                oSetUniqueCoords.end())
            {
                oSetUniqueCoords.insert(std::pair<int, int>(nX, nY));

                // Calculate the difference from the last point and add encoded coordinates
                int nDiffX = nX - nLastX;
                int nDiffY = nY - nLastY;
                poGPBFeature->addGeometry(EncodeSInt(nDiffX));
                poGPBFeature->addGeometry(EncodeSInt(nDiffY));

                // Update last point coordinates
                nLastX = nX;
                nLastY = nY;
            }
        }
    }

    // If we have unique points, update the geometry and return success
    GUInt32 nPoints = static_cast<GUInt32>(oSetUniqueCoords.size());
    bool bGeomOK = nPoints > 0;
    poGPBFeature->setGeometry(0, GetCmdCountCombined(knCMD_MOVETO, nPoints));
    return bGeomOK;
}

/************************************************************************/
/*                          HandleLineStringGeometry()                  */
/************************************************************************/
bool GeometryEncoder::HandleLineStringGeometry(
    const OGRLineString *poLS, MVTTileLayerFeature *poGPBFeature,
    double &dfAreaOrLength) const
{
    std::cout << "OGRMVTWriterDataset::HandleLineStringGeometry" << std::endl;

    // Constants specific to this method
    const bool bWriteLastPoint = true;
    const bool bReverseOrder = false;
    const GUInt32 nMinLineTo = 1;

    int nLastX = 0;
    int nLastY = 0;
    OGRLineString oOutLS;

    // Call the EncodeLineString function
    bool bGeomOK =
        EncodeLineString(poGPBFeature, poLS, &oOutLS, bWriteLastPoint,
                         bReverseOrder, nMinLineTo, nLastX, nLastY);

    // Update the length of the line string
    if (bGeomOK)
    {
        dfAreaOrLength = oOutLS.get_Length();
    }

    return bGeomOK;
}

/************************************************************************/
/*                          HandleMultiLineStringOrCollection()         */
/************************************************************************/
bool GeometryEncoder::HandleMultiLineStringOrCollection(
    const OGRGeometryCollection *poGC, MVTTileLayerFeature *poGPBFeature,
    double &dfAreaOrLength) const
{
    std::cout << "HandleMultiLineStringOrCollection" << std::endl;

    // Constants specific to this method
    const bool bWriteLastPoint = true;
    const bool bReverseOrder = false;
    const GUInt32 nMinLineTo = 1;

    bool bGeomOK = false;
    int nLastX = 0;
    int nLastY = 0;

    // Iterate through the sub-geometries in the collection
    for (auto &&poSubGeom : poGC)
    {
        if (wkbFlatten(poSubGeom->getGeometryType()) == wkbLineString)
        {
            const OGRLineString *poLS = poSubGeom->toLineString();
            OGRLineString oOutLS;

            // Encode each line string in the collection
            bool bSubGeomOK =
                EncodeLineString(poGPBFeature, poLS, &oOutLS, bWriteLastPoint,
                                 bReverseOrder, nMinLineTo, nLastX, nLastY);

            if (bSubGeomOK)
            {
                dfAreaOrLength +=
                    oOutLS.get_Length();  // Accumulate total length
            }

            bGeomOK |= bSubGeomOK;  // Update geometry validity
        }
    }

    return bGeomOK;
}

/************************************************************************/
/*                          HandlePolygonGeometry()                     */
/************************************************************************/
bool GeometryEncoder::HandlePolygonGeometry(const OGRPolygon *poPoly,
                                            MVTTileLayerFeature *poGPBFeature,
                                            double &dfArea) const
{
    std::cout << "HandlePolygonGeometry" << std::endl;

    int nLastX = 0;
    int nLastY = 0;
    OGRPolygon oOutPoly;
    const GUInt32 nInitialSize = poGPBFeature->getGeometryCount();

    // Encode the polygon
    bool bGeomOK =
        EncodePolygon(poGPBFeature, poPoly, &oOutPoly, nLastX, nLastY, dfArea);

    // Check validity of the resulting polygon
    int bIsValid;
    {
        CPLErrorStateBackuper oErrorStateBackuper(CPLQuietErrorHandler);
        bIsValid = oOutPoly.IsValid();
    }

    if (!bIsValid)
    {
        // Make the polygon valid and emit it
        std::unique_ptr<OGRGeometry> poPolyValid(oOutPoly.MakeValid());
        if (poPolyValid)
        {
            poGPBFeature->resizeGeometryArray(nInitialSize);
            EmitValidPolygon(poPolyValid.get(), poGPBFeature, dfArea);
        }
    }

    return bGeomOK;
}

/************************************************************************/
/*                          HandleMultiPolygonOrCollection()            */
/************************************************************************/
bool GeometryEncoder::HandleMultiPolygonOrCollection(
    const OGRGeometryCollection *poGC, MVTTileLayerFeature *poGPBFeature,
    double &dfTotalArea) const
{
    std::cout << "HandleMultiPolygonOrCollection" << std::endl;

    int nLastX = 0;
    int nLastY = 0;
    OGRMultiPolygon oOutMP;
    const GUInt32 nInitialSize = poGPBFeature->getGeometryCount();
    bool bGeomOK = false;

    for (auto &&poSubGeom : poGC)
    {
        if (wkbFlatten(poSubGeom->getGeometryType()) == wkbPolygon)
        {
            const OGRPolygon *poPoly = poSubGeom->toPolygon();
            double dfPartArea = 0.0;
            auto poOutPoly = std::make_unique<OGRPolygon>();

            // Encode each polygon in the collection
            bool bSubGeomOK =
                EncodePolygon(poGPBFeature, poPoly, poOutPoly.get(), nLastX,
                              nLastY, dfPartArea);

            if (bSubGeomOK)
            {
                dfTotalArea += dfPartArea;
                oOutMP.addGeometryDirectly(poOutPoly.release());
            }

            bGeomOK |= bSubGeomOK;
        }
    }

    // Check validity of the resulting multipolygon
    int bIsValid;
    {
        CPLErrorStateBackuper oErrorStateBackuper(CPLQuietErrorHandler);
        bIsValid = oOutMP.IsValid();
    }

    if (!bIsValid)
    {
        // Make the multipolygon valid and emit it
        std::unique_ptr<OGRGeometry> poMPValid(oOutMP.MakeValid());
        if (poMPValid)
        {
            poGPBFeature->resizeGeometryArray(nInitialSize);
            EmitValidPolygon(poMPValid.get(), poGPBFeature, dfTotalArea);
        }
    }

    return bGeomOK;
}

/************************************************************************/
/*                          EncodeLineString()                          */
/************************************************************************/
bool GeometryEncoder::EncodeLineString(MVTTileLayerFeature *poGPBFeature,
                                       const OGRLineString *poLS,
                                       OGRLineString *poOutLS,
                                       bool bWriteLastPoint, bool bReverseOrder,
                                       GUInt32 nMinLineTo, int &nLastX,
                                       int &nLastY, double dfTopX,
                                       double dfTopY, double dfTileDim) const
{
    std::cout << "EncodeLineString 1" << std::endl;

    const GUInt32 nInitialSize = poGPBFeature->getGeometryCount();
    const int nLastXOri = nLastX;
    const int nLastYOri = nLastY;
    GUInt32 nLineToCount = 0;
    const int nPoints = poLS->getNumPoints() - (bWriteLastPoint ? 0 : 1);
    if (poOutLS)
        poOutLS->setNumPoints(nPoints);
    int nFirstX = 0;
    int nFirstY = 0;
    int nLastXValid = nLastX;
    int nLastYValid = nLastY;
    for (int i = 0; i < nPoints; i++)
    {
        int nX, nY;
        int nSrcIdx = bReverseOrder ? poLS->getNumPoints() - 1 - i : i;
        double dfX = poLS->getX(nSrcIdx);
        double dfY = poLS->getY(nSrcIdx);
        ConvertToTileCoords(dfX, dfY, nX, nY, dfTopX, dfTopY, dfTileDim);
        int nDiffX = nX - nLastX;
        int nDiffY = nY - nLastY;
        if (i == 0 || nDiffX != 0 || nDiffY != 0)
        {
            if (i > 0)
            {
                nLineToCount++;
                if (nLineToCount == 1)
                {
                    poGPBFeature->addGeometry(
                        GetCmdCountCombined(knCMD_MOVETO, 1));
                    const int nLastDiffX = nLastX - nLastXOri;
                    const int nLastDiffY = nLastY - nLastYOri;
                    poGPBFeature->addGeometry(EncodeSInt(nLastDiffX));
                    poGPBFeature->addGeometry(EncodeSInt(nLastDiffY));
                    if (poOutLS)
                        poOutLS->setPoint(0, nLastX, nLastY);

                    // To be modified later
                    poGPBFeature->addGeometry(
                        GetCmdCountCombined(knCMD_LINETO, 0));
                }

                poGPBFeature->addGeometry(EncodeSInt(nDiffX));
                poGPBFeature->addGeometry(EncodeSInt(nDiffY));
                if (poOutLS)
                    poOutLS->setPoint(nLineToCount, nX, nY);
            }
            else
            {
                nFirstX = nX;
                nFirstY = nY;
            }
            nLastXValid = nLastX;
            nLastYValid = nLastY;
            nLastX = nX;
            nLastY = nY;
        }
    }

    // If last point of ring is identical to first one, discard it
    if (nMinLineTo == 2 && nLineToCount > 0 && nFirstX == nLastX &&
        nFirstY == nLastY)
    {
        poGPBFeature->resizeGeometryArray(poGPBFeature->getGeometryCount() - 2);
        nLineToCount--;
        nLastX = nLastXValid;
        nLastY = nLastYValid;
    }

    if (nLineToCount >= nMinLineTo)
    {
        if (poOutLS)
            poOutLS->setNumPoints(1 + nLineToCount);
        // Patch actual number of points in LINETO command
        poGPBFeature->setGeometry(
            nInitialSize + 3, GetCmdCountCombined(knCMD_LINETO, nLineToCount));
        return true;
    }
    else
    {
        poGPBFeature->resizeGeometryArray(nInitialSize);
        nLastX = nLastXOri;
        nLastY = nLastYOri;
        return false;
    }
}

/************************************************************************/
/*                          EncodeLineString()                          */
/************************************************************************/
bool GeometryEncoder::EncodeLineString(MVTTileLayerFeature *poGPBFeature,
                                       const OGRLineString *poLS,
                                       OGRLineString *poOutLS,
                                       bool bWriteLastPoint, bool bReverseOrder,
                                       GUInt32 nMinLineTo, int &nLastX,
                                       int &nLastY) const
{
    std::cout << "EncodeLineString 2" << std::endl;

    return EncodeLineString(poGPBFeature, poLS, poOutLS, bWriteLastPoint,
                            bReverseOrder, nMinLineTo, nLastX, nLastY, m_dfTopX,
                            m_dfTopY, m_dfTileDim);
}

/************************************************************************/
/*                          EncodePolygon()                             */
/************************************************************************/
bool GeometryEncoder::EncodePolygon(MVTTileLayerFeature *poGPBFeature,
                                    const OGRPolygon *poPoly,
                                    OGRPolygon *poOutPoly, int &nLastX,
                                    int &nLastY, double &dfArea, double dfTopX,
                                    double dfTopY, double dfTileDim) const
{
    std::cout << "EncodePolygon 1" << std::endl;

    dfArea = 0;
    auto poOutOuterRing = std::make_unique<OGRLinearRing>();
    for (int i = 0; i < 1 + poPoly->getNumInteriorRings(); i++)
    {
        const OGRLinearRing *poRing = (i == 0) ? poPoly->getExteriorRing()
                                               : poPoly->getInteriorRing(i - 1);
        if (poRing->getNumPoints() < 4 ||
            poRing->getX(0) != poRing->getX(poRing->getNumPoints() - 1) ||
            poRing->getY(0) != poRing->getY(poRing->getNumPoints() - 1))
        {
            if (i == 0)
                return false;
            continue;
        }
        const bool bWriteLastPoint = false;
        // If dealing with input geometry in CRS units, exterior rings must
        // be clockwise oriented.
        // But if re-encoding a geometry already in tile coordinates
        // (dfTileDim == 0), this is the reverse.
        const bool bReverseOrder = dfTileDim != 0
                                       ? ((i == 0 && !poRing->isClockwise()) ||
                                          (i > 0 && poRing->isClockwise()))
                                       : ((i == 0 && poRing->isClockwise()) ||
                                          (i > 0 && !poRing->isClockwise()));
        const GUInt32 nMinLineTo = 2;
        std::unique_ptr<OGRLinearRing> poOutInnerRing;
        if (i > 0)
            poOutInnerRing = std::make_unique<OGRLinearRing>();
        OGRLinearRing *poOutRing =
            poOutInnerRing.get() ? poOutInnerRing.get() : poOutOuterRing.get();

        bool bSuccess = EncodeLineString(
            poGPBFeature, poRing, poOutRing, bWriteLastPoint, bReverseOrder,
            nMinLineTo, nLastX, nLastY, dfTopX, dfTopY, dfTileDim);
        if (!bSuccess)
        {
            if (i == 0)
                return false;
            continue;
        }

        if (poOutPoly == nullptr)
        {
            poGPBFeature->addGeometry(GetCmdCountCombined(knCMD_CLOSEPATH, 1));
            continue;
        }

        poOutRing->closeRings();

        poOutPoly->addRing(poOutRing);
        if (i > 0)
            dfArea -= poOutRing->get_Area();
        else
            dfArea = poOutRing->get_Area();

        poGPBFeature->addGeometry(GetCmdCountCombined(knCMD_CLOSEPATH, 1));
    }

    return true;
}

/************************************************************************/
/*                          EncodePolygon()                             */
/************************************************************************/
bool GeometryEncoder::EncodePolygon(MVTTileLayerFeature *poGPBFeature,
                                    const OGRPolygon *poPoly,
                                    OGRPolygon *poOutPoly, int &nLastX,
                                    int &nLastY, double &dfArea) const
{
    std::cout << "EncodePolygon 2" << std::endl;

    return EncodePolygon(poGPBFeature, poPoly, poOutPoly, nLastX, nLastY,
                         dfArea, m_dfTopX, m_dfTopY, m_dfTileDim);
}

/************************************************************************/
/*                     EncodeRepairedOuterRing()                        */
/************************************************************************/
bool GeometryEncoder::EncodeRepairedOuterRing(MVTTileLayerFeature *poGPBFeature,
                                              const OGRPolygon &oInPoly,
                                              int &nLastX, int &nLastY) const
{
    std::cout << "EncodeRepairedOuterRing" << std::endl;

    std::unique_ptr<OGRGeometry> poFixedGeom(oInPoly.Buffer(0));
    if (!poFixedGeom.get() || poFixedGeom->IsEmpty())
    {
        return false;
    }

    OGRPolygon *poPoly = nullptr;
    if (wkbFlatten(poFixedGeom->getGeometryType()) == wkbMultiPolygon)
    {
        OGRMultiPolygon *poMP = poFixedGeom.get()->toMultiPolygon();
        poPoly = poMP->getGeometryRef(0)->toPolygon();
    }
    else if (wkbFlatten(poFixedGeom->getGeometryType()) == wkbPolygon)
    {
        poPoly = poFixedGeom.get()->toPolygon();
    }
    if (!poPoly)
        return false;

    OGRLinearRing *poRing = poPoly->getExteriorRing();
    const bool bReverseOrder = !poRing->isClockwise();

    const GUInt32 nInitialSize = poGPBFeature->getGeometryCount();
    const int nLastXOri = nLastX;
    const int nLastYOri = nLastY;
    GUInt32 nLineToCount = 0;
    const int nPoints = poRing->getNumPoints() - 1;
    auto poOutLinearRing = std::make_unique<OGRLinearRing>();
    poOutLinearRing->setNumPoints(nPoints);
    for (int i = 0; i < nPoints; i++)
    {
        int nSrcIdx = bReverseOrder ? poRing->getNumPoints() - 1 - i : i;
        double dfX = poRing->getX(nSrcIdx);
        double dfY = poRing->getY(nSrcIdx);
        int nX = static_cast<int>(std::round(dfX));
        int nY = static_cast<int>(std::round(dfY));
        if (nX != dfX || nY != dfY)
            continue;
        int nDiffX = nX - nLastX;
        int nDiffY = nY - nLastY;
        if (i == 0 || nDiffX != 0 || nDiffY != 0)
        {
            if (i > 0)
            {
                nLineToCount++;
                if (nLineToCount == 1)
                {
                    poGPBFeature->addGeometry(
                        GetCmdCountCombined(knCMD_MOVETO, 1));
                    const int nLastDiffX = nLastX - nLastXOri;
                    const int nLastDiffY = nLastY - nLastYOri;
                    poGPBFeature->addGeometry(EncodeSInt(nLastDiffX));
                    poGPBFeature->addGeometry(EncodeSInt(nLastDiffY));
                    poOutLinearRing->setPoint(0, nLastX, nLastY);

                    // To be modified later
                    poGPBFeature->addGeometry(
                        GetCmdCountCombined(knCMD_LINETO, 0));
                }

                poGPBFeature->addGeometry(EncodeSInt(nDiffX));
                poGPBFeature->addGeometry(EncodeSInt(nDiffY));
                poOutLinearRing->setPoint(nLineToCount, nX, nY);
            }
            nLastX = nX;
            nLastY = nY;
        }
    }
    if (nLineToCount >= 2)
    {
        poOutLinearRing->setNumPoints(1 + nLineToCount);
        OGRPolygon oOutPoly;
        oOutPoly.addRingDirectly(poOutLinearRing.release());
        int bIsValid;
        {
            CPLErrorStateBackuper oErrorStateBackuper(CPLQuietErrorHandler);
            bIsValid = oOutPoly.IsValid();
        }
        if (bIsValid)
        {
            // Patch actual number of points in LINETO command
            poGPBFeature->setGeometry(
                nInitialSize + 3,
                GetCmdCountCombined(knCMD_LINETO, nLineToCount));
            poGPBFeature->addGeometry(GetCmdCountCombined(knCMD_CLOSEPATH, 1));
            return true;
        }
    }

    poGPBFeature->resizeGeometryArray(nInitialSize);
    nLastX = nLastXOri;
    nLastY = nLastYOri;
    return false;
}

/************************************************************************/
/*                          EmitValidPolygon()                          */
/************************************************************************/
bool GeometryEncoder::EmitValidPolygon(const OGRGeometry *poValidGeom,
                                       MVTTileLayerFeature *poGPBFeature,
                                       double &dfAreaOrLength) const
{
    std::cout << "EmitValidPolygon" << std::endl;

    bool bGeomOK = false;
    dfAreaOrLength = 0.0;
    int nLastX = 0;
    int nLastY = 0;

    if (wkbFlatten(poValidGeom->getGeometryType()) == wkbPolygon)
    {
        const OGRPolygon *poPoly = poValidGeom->toPolygon();
        double dfPartArea = 0.0;
        bGeomOK = EncodePolygon(poGPBFeature, poPoly, nullptr, nLastX, nLastY,
                                dfPartArea, 0, 0, 0);
        dfAreaOrLength = dfPartArea;
    }
    else if (OGR_GT_IsSubClassOf(poValidGeom->getGeometryType(),
                                 wkbGeometryCollection))
    {
        for (auto &&poSubGeom : poValidGeom->toGeometryCollection())
        {
            if (wkbFlatten(poSubGeom->getGeometryType()) == wkbPolygon)
            {
                const OGRPolygon *poPoly = poSubGeom->toPolygon();
                double dfPartArea = 0.0;
                bGeomOK |= EncodePolygon(poGPBFeature, poPoly, nullptr, nLastX,
                                         nLastY, dfPartArea, 0, 0, 0);
                dfAreaOrLength += dfPartArea;
            }
            else if (wkbFlatten(poSubGeom->getGeometryType()) ==
                     wkbMultiPolygon)
            {
                const OGRMultiPolygon *poMPoly = poSubGeom->toMultiPolygon();
                for (const auto *poPoly : poMPoly)
                {
                    double dfPartArea = 0.0;
                    bGeomOK |=
                        EncodePolygon(poGPBFeature, poPoly, nullptr, nLastX,
                                      nLastY, dfPartArea, 0, 0, 0);
                    dfAreaOrLength += dfPartArea;
                }
            }
        }
    }

    return bGeomOK;
}
