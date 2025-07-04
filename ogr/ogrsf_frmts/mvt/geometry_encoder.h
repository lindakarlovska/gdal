/******************************************************************************
 *
 * Project:  MVT Translator
 * Purpose:  Geometry clipping, transformation and encoding to MVT tile features.
 * Author:   Even Rouault, Even Rouault <even dot rouault at spatialys dot com>
 *           Linda Karlovska <linda dot karlovska at seznam dot cz>
 *
 ******************************************************************************
 * Copyright (c) 2018, Even Rouault <even dot rouault at spatialys dot com>
 * Copyright (c) 2025, Linda Karlovska <linda dot karlovska at seznam dot cz>
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#ifndef GEOMETRY_ENCODER_H
#define GEOMETRY_ENCODER_H

#include "ogr_geometry.h"
#include "mvt_tile.h"
#include "gpb.h"

/**
 * @class GeometryEncoder
 * @brief Encodes OGR geometries into Mapbox Vector Tile (MVT) protobuf features.
 *
 * The encoded features are stored in a temporary SQLite database.
 * The temporary database accumulates features for multiple tiles.
 * Later, entire tiles are assembled from the temporary storage
 * and encoded into final PBF (protobuf) files or MBTiles database.
 *
 * Note:
 * Compression and final storage are handled outside this class.
 */

class GeometryEncoder
{
  public:
    /**
     * @brief Constructor.
     * 
     * @param dfTileDim   Dimension (width/height) of the tile in projected coordinate units.
     * @param dfTopX      X coordinate of the top-left corner of the tile (projected space).
     * @param dfTopY      Y coordinate of the top-left corner of the tile (projected space).
     * @param nExtent     MVT tile extent (typically 4096).
     */
    GeometryEncoder(double dfTileDim, double dfTopX, double dfTopY,
                    int nExtent);

    /**
     * @brief Encodes a geometry into MVT feature format.
     * @param poOriginalGeom Original geometry for fallback emission.
     * @param poGeomToEncode Clipped and possibly simplified geometry.
     * @param poGPBFeature Output MVT feature.
     * @param dfAreaOrLength Output area (for polygons) or length (for lines).
     * @return True if encoding succeeded.
     */

    bool Encode(const OGRGeometry *poOriginalGeom,
                const OGRGeometry *poGeomToEncode,
                MVTTileLayerFeature *poGPBFeature,
                double &dfAreaOrLength) const;

  private:
    double m_dfTileDim;
    double m_dfTopX;
    double m_dfTopY;
    int m_nExtent;

	/**
	 * @brief Converts map coordinates to tile-relative coordinates.
	 *
	 * This version allows explicit control over the tile origin and size.
	 * It converts real-world coordinates (e.g., in meters or degrees)
	 * into local tile coordinates, using a given tile top-left origin and tile dimension.
	 *
	 * @param dfX        World X coordinate (e.g., WebMercator meters)
	 * @param dfY        World Y coordinate
	 * @param[out] nX    Output X in tile coordinates (integer pixel)
	 * @param[out] nY    Output Y in tile coordinates
	 * @param dfTopX     Top-left X coordinate of the tile
	 * @param dfTopY     Top-left Y coordinate of the tile
	 * @param dfTileDim  Tile dimension in world units (e.g., width/height in meters)
	 */
    void ConvertToTileCoords(double dfX, double dfY, int &nX, int &nY,
                             double dfTopX, double dfTopY,
                             double dfTileDim) const;

    /**
	 * @brief Converts map coordinates to tile-relative coordinates using internal tile parameters.
	 *
	 * This overload uses the tile origin and tile dimension defined internally in the encoder.
	 * Converts real-world coordinates (e.g., in meters or degrees) into tile-local coordinates.
	 *
	 * @param dfX     World X coordinate (e.g., WebMercator meters)
	 * @param dfY     World Y coordinate
	 * @param[out] nX Output X in tile coordinates (integer pixel)
	 * @param[out] nY Output Y in tile coordinates
	 */
    void ConvertToTileCoords(double dfX, double dfY, int &nX, int &nY) const;

    /**
	 * @brief Handles encoding of a single point geometry.
	 * 
	 * @param poPoint The input point.
	 * @param poGPBFeature Output PBF feature.
	 * @param oTileParams Tile parameters.
	 * @return true if successful.
	 */
    bool HandlePointGeometry(const OGRPoint *poPoint,
                             MVTTileLayerFeature *poGPBFeature) const;

    /**
	 * @brief Handles encoding of a multipoint collection.
	 * 
	 * @param poGC The input geometry collection.
	 * @param poGPBFeature Output PBF feature.
	 * @param oTileParams Tile parameters.
	 * @return true if successful.
	 */
    bool HandleMultiPointOrCollection(const OGRGeometryCollection *poGC,
                                      MVTTileLayerFeature *poGPBFeature) const;

    /**
	 * @brief Handles encoding of a LineString geometry.
	 * 
	 * @param poLS The input LineString.
	 * @param poGPBFeature Output PBF feature.
	 * @param oTileParams Tile parameters.
	 * @param dfAreaOrLength Output length of the line.
	 * @return true if successful.
	 */
    bool HandleLineStringGeometry(const OGRLineString *poLS,
                                  MVTTileLayerFeature *poGPBFeature,
                                  double &dfAreaOrLength) const;

    /**
	 * @brief Handles encoding of a MultiLineString or geometry collection.
	 * 
	 * @param poGC The input geometry collection.
	 * @param poGPBFeature Output PBF feature.
	 * @param oTileParams Tile parameters.
	 * @param dfAreaOrLength Output total length.
	 * @return true if successful.
	 */
    bool HandleMultiLineStringOrCollection(const OGRGeometryCollection *poGC,
                                           MVTTileLayerFeature *poGPBFeature,
                                           double &dfAreaOrLength) const;

    /**
	 * @brief Handles encoding of a Polygon geometry.
	 * 
	 * @param poPoly The input polygon.
	 * @param poGPBFeature Output PBF feature.
	 * @param oTileParams Tile parameters.
	 * @param dfArea Output area of the polygon.
	 * @return true if successful.
	 */
    bool HandlePolygonGeometry(const OGRPolygon *poPoly,
                               MVTTileLayerFeature *poGPBFeature,
                               double &dfArea) const;

    /**
	 * @brief Handles encoding of a MultiPolygon or geometry collection.
	 * 
	 * @param poGC The input geometry collection.
	 * @param poGPBFeature Output PBF feature.
	 * @param oTileParams Tile parameters.
	 * @param dfTotalArea Output total area of the multipolygon.
	 * @return true if successful.
	 */
    bool HandleMultiPolygonOrCollection(const OGRGeometryCollection *poGC,
                                        MVTTileLayerFeature *poGPBFeature,
                                        double &dfTotalArea) const;

    /**
     * @brief Encodes a LineString into the MVT feature format.
     *
     * This variant uses explicit tile transformation parameters to convert
     * the geometry coordinates into local tile space.
     *
     * @param poGPBFeature Output MVT feature to store the encoded geometry.
     * @param poLS Input LineString geometry to encode.
     * @param poOutLS Optional pointer to an output LineString receiving the transformed geometry.
     *                Can be nullptr if not needed.
     * @param bWriteLastPoint Whether to explicitly include the last point in the output.
     * @param bReverseOrder Whether to reverse the point order before encoding.
     * @param nMinLineTo Minimum number of points required to emit a LineTo command.
     * @param nLastX Reference to last encoded X tile coordinate (used for delta encoding).
     * @param nLastY Reference to last encoded Y tile coordinate.
     * @param dfTopX Left/top tile X coordinate (tile origin).
     * @param dfTopY Left/top tile Y coordinate (tile origin).
     * @param dfTileDim Size of the tile in coordinate units (extent).
     * @return true if the encoding was successful and the geometry met the required conditions.
     */
    bool EncodeLineString(MVTTileLayerFeature *poGPBFeature,
                          const OGRLineString *poLS, OGRLineString *poOutLS,
                          bool bWriteLastPoint, bool bReverseOrder,
                          GUInt32 nMinLineTo, int &nLastX, int &nLastY,
                          double dfTopX, double dfTopY, double dfTileDim) const;

    /**
     * @brief Encodes a LineString using internal tile transformation parameters.
     *
     * Convenience overload that uses the tile parameters stored in the
     * GeometryEncoder instance. Appropriate when default transformation is desired.
     *
     * @param poGPBFeature Output MVT feature to store the encoded geometry.
     * @param poLS Input LineString geometry to encode.
     * @param poOutLS Optional pointer to an output LineString receiving the transformed geometry.
     * @param bWriteLastPoint Whether to explicitly include the last point in the output.
     * @param bReverseOrder Whether to reverse the point order before encoding.
     * @param nMinLineTo Minimum number of points required to emit a LineTo command.
     * @param nLastX Reference to last encoded X tile coordinate (used for delta encoding).
     * @param nLastY Reference to last encoded Y tile coordinate.
     */
    bool EncodeLineString(MVTTileLayerFeature *poGPBFeature,
                          const OGRLineString *poLS, OGRLineString *poOutLS,
                          bool bWriteLastPoint, bool bReverseOrder,
                          GUInt32 nMinLineTo, int &nLastX, int &nLastY) const;

    /**
	 * @brief Encodes a polygon into the MVT feature format.
	 * 
	 * This variant allows passing explicit tile transformation parameters, which are
         * used to convert the input polygon coordinates into tile-local coordinates.
         *
     * @param poGPBFeature Output MVT feature to which the geometry will be encoded.
     * @param poPoly Input polygon geometry to encode.
     * @param poOutPoly Output polygon that receives the transformed geometry. Can be nullptr.
     * @param nLastX Reference to last X tile coordinate; updated during encoding.
     * @param nLastY Reference to last Y tile coordinate; updated during encoding.
     * @param dfArea Reference to a variable receiving the computed area of the polygon.
     * @param dfTopX Left/top tile X coordinate (tile origin).
     * @param dfTopY Left/top tile Y coordinate (tile origin).
     * @param dfTileDim Size of the tile in coordinate units.
     * @return true if the encoding was successful and the geometry was valid.
	 */
    bool EncodePolygon(MVTTileLayerFeature *poGPBFeature,
                       const OGRPolygon *poPoly, OGRPolygon *poOutPoly,
                       int &nLastX, int &nLastY, double &dfArea, double dfTopX,
                       double dfTopY, double dfTileDim) const;

    /**
     * @brief Encodes a polygon using internal tile transformation parameters.
     *
     * This is a convenience overload that uses the tile parameters stored in the
     * GeometryEncoder instance. Useful when no custom transformation is required.
     *
     * @param poGPBFeature Output MVT feature to which the geometry will be encoded.
     * @param poPoly Input polygon geometry to encode.
     * @param poOutPoly Output polygon that receives the transformed geometry. Can be nullptr.
     * @param nLastX Reference to last X tile coordinate; updated during encoding.
     * @param nLastY Reference to last Y tile coordinate; updated during encoding.
     * @param dfArea Reference to a variable receiving the computed area of the polygon.
     * @return true if the encoding was successful and the geometry was valid.
     */
    bool EncodePolygon(MVTTileLayerFeature *poGPBFeature,
                       const OGRPolygon *poPoly, OGRPolygon *poOutPoly,
                       int &nLastX, int &nLastY, double &dfArea) const;

    /**
	 * @brief Attempts to encode a repaired version of the outer ring of a polygon into a MVT feature.
	 *
	 * This function performs a repair on the polygon using a 0-buffer (which fixes self-intersections and other invalidities).
	 * It extracts the outer ring, ensures it is oriented clockwise, and encodes it as a series of MVT commands 
	 * (MOVE_TO, LINE_TO, CLOSE_PATH) into the given tile feature.

	 * @param poGPBFeature Output PBF feature to which the geometry will be added.
	 * @param oInPoly Input polygon to repair and encode.
	 * @param nLastX Last used X coordinate in tile space (input/output, used for delta encoding).
	 * @param nLastY Last used Y coordinate in tile space (input/output, used for delta encoding).
	 * @return true if a valid and non-degenerate outer ring was encoded successfully; false otherwise.
	 */
    bool EncodeRepairedOuterRing(MVTTileLayerFeature *poGPBFeature,
                                 const OGRPolygon &oInPoly, int &nLastX,
                                 int &nLastY) const;

    /**
	 * @brief Emits a valid polygon geometry into an MVT tile feature.
	 * 
	 * @param poValidGeom       Pointer to a valid OGRGeometry (Polygon or GeometryCollection)
	 * @param poGPBFeature      Pointer to the target MVT feature object for encoding
	 * @param dfAreaOrLength    Output reference, will contain the total area of encoded polygons
	 *
	 * @return true if at least one polygon was successfully encoded and added to the feature, false otherwise
	 *
	 */
    bool EmitValidPolygon(const OGRGeometry *poValidGeom,
                          MVTTileLayerFeature *poGPBFeature,
                          double &dfAreaOrLength) const;
};

#endif  // GEOMETRY_ENCODER_H
