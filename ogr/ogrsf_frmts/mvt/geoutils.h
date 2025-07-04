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

// WebMercator related constants
constexpr double kmSPHERICAL_RADIUS = 6378137.0;

/**
 * @brief Utility class for coordinate and envelope transformations related to MVT.
 */
class MVTGeoUtils
{
  public:
    /**
    * @brief Initializes the spatial reference and tiling parameters for Web Mercator.
    *
    * @param poSRS Pointer to an OGRSpatialReference object to be initialized to EPSG:3857.
    * @param[out] dfTopX Reference to a double that will receive the X coordinate of the top-left corner.
    * @param[out] dfTopY Reference to a double that will receive the Y coordinate of the top-left corner.
    * @param[out] dfTileDim0 Reference to a double that will receive the tile dimension (extent size).
    */
    static void InitWebMercatorTilingScheme(OGRSpatialReference *poSRS,
                                            double &dfTopX, double &dfTopY,
                                            double &dfTileDim0);
    /**
     * @brief Converts coordinates from Spherical Mercator projection to longitude and latitude (WGS84).
     *
     * @param[in,out] pdfX Pointer to the X coordinate (meters in Spherical Mercator).
     * @param[in,out] pdfY Pointer to the Y coordinate (meters in Spherical Mercator).
     */
    static void SphericalMercatorToLongLat(double *pdfX, double *pdfY);

    /**
     * @brief Converts longitude and latitude (WGS84) to coordinates in Spherical Mercator projection.
     *
     * @param[in,out] pdfX Pointer to longitude in degrees.
     * @param[in,out] pdfY Pointer to latitude in degrees.
     */
    static void LongLatToSphericalMercator(double *pdfX, double *pdfY);

    /**
     * @brief Converts an envelope from WGS84 coordinates to the target spatial reference system (SRS).
     *
     * If the target SRS is EPSG:3857 (Spherical Mercator), a direct conversion is performed,
     * otherwise, an OGRCoordinateTransformation is used.
     *
     * @param[in] poTargetSRS Target spatial reference system.
     * @param[in,out] dfX0 Left lower corner X coordinate.
     * @param[in,out] dfY0 Left lower corner Y coordinate.
     * @param[in,out] dfX1 Right upper corner X coordinate.
     * @param[in,out] dfY1 Right upper corner Y coordinate.
     */
    static void ConvertFromWGS84(OGRSpatialReference *poTargetSRS, double &dfX0,
                                 double &dfY0, double &dfX1, double &dfY1);

    /**
     * @brief Converts an envelope from Spherical Mercator projection to WGS84 longitude and latitude,
     * and clamps the latitude values to the range [-85, 85].
     *
     * @param[in] oEnvelope Input envelope in Spherical Mercator.
     * @param[out] oTransformedEnvelope Output envelope in WGS84 coordinates.
     */
    static void ComputeStandardEnvelope(const OGREnvelope &oEnvelope,
                                        OGREnvelope &oTransformedEnvelope);

    /**
     * @brief Transforms an envelope from a given spatial reference system to WGS84 (EPSG:4326).
     *
     * Transforms all four corners of the envelope and computes the bounding envelope of the transformed points.
     *
     * @param[in] oEnvelope Input envelope.
     * @param[in] poSRS Spatial reference system of the input envelope.
     * @param[out] oTransformedEnvelope Output envelope in WGS84 coordinates.
     */
    static void TransformEnvelopeToWGS84(const OGREnvelope &oEnvelope,
                                         OGRSpatialReference *poSRS,
                                         OGREnvelope &oTransformedEnvelope);

    /**
     * @brief Computes the center of the envelope and returns it as a string in the format "longitude,latitude,minZoom".
     *
     * @param[in] oEnvelope Input envelope.
     * @param[in] nMinZoom Minimum zoom level.
     * @return CPLString String representing the center and minimum zoom.
     */
    static CPLString ComputeCenter(const OGREnvelope &oEnvelope, int nMinZoom);

    /**
     * @brief Computes the bounding box of the envelope and returns it as a string in the format "minX,minY,maxX,maxY".
     *
     * @param[in] oEnvelope Input envelope.
     * @return CPLString String representing the bounding box.
     */
    static CPLString ComputeBounds(const OGREnvelope &oEnvelope);
};

#endif  // GEOUTILS_H
