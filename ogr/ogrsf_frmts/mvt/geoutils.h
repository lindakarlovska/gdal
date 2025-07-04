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

#ifndef GEOUTILS_H
#define GEOUTILS_H

#include "cpl_string.h"
#include "ogrsf_frmts.h"

/**
 * @brief Utility namespace for coordinate and envelope transformations related to MVT.
 */
namespace MVTGeoUtils
{
// Constants for WebMercator and geographic (WGS84) coordinate systems
constexpr double kmSPHERICAL_RADIUS = 6378137.0;
constexpr double kmMAX_GM = kmSPHERICAL_RADIUS * M_PI;
constexpr const char *kSRSWebMercator = SRS_WKT_WGS84_PSEUDO_MERCATOR;
constexpr const char *kSRSWGS84 = SRS_WKT_WGS84_LAT_LONG;

void InitWebMercatorTilingScheme(OGRSpatialReference *poSRS, double &dfTopX,
                                 double &dfTopY, double &dfTileDim0);

bool InitUserDefinedTilingScheme(const char *pszTilingScheme,
                                 OGRSpatialReference *poSRS, double &dfTopX,
                                 double &dfTopY, double &dfTileDim0,
                                 int &nTileMatrixWidth0,
                                 int &nTileMatrixHeight0);

void SphericalMercatorToLongLat(double *pdfX, double *pdfY);
void LongLatToSphericalMercator(double *pdfX, double *pdfY);
void ConvertFromWGS84(OGRSpatialReference *poTargetSRS, double &dfX0,
                      double &dfY0, double &dfX1, double &dfY1);

void ComputeStandardEnvelope(const OGREnvelope &oEnvelope,
                             OGREnvelope &oTransformedEnvelope);
void TransformEnvelopeToWGS84(const OGREnvelope &oEnvelope,
                              OGRSpatialReference *poSRS,
                              OGREnvelope &oTransformedEnvelope);

CPLString ComputeCenter(const OGREnvelope &oEnvelope, int nMinZoom);
CPLString ComputeBounds(const OGREnvelope &oEnvelope);

bool ValidateMinMaxZoom(int nMinZoom, int nMaxZoom);
};  // namespace MVTGeoUtils

#endif  // GEOUTILS_H
