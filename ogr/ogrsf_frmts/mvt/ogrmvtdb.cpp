/******************************************************************************
 *
 * Project:  MVT Translator
 * Purpose:  Mapbox Vector Tile (MVT) decoder and encoder
 * Authors:  Even Rouault, <even dot rouault at spatialys dot com>
 *           Linda Karlovska, <linda dot karlovska at seznam dot cz> (refactoring)
 *
 ******************************************************************************
 * Copyright (c) 2018, Even Rouault <even dot rouault at spatialys dot com>
 *               2025, Linda Karlovska <linda dot karlovska at seznam dot cz>
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#include "ogrmvtdb.h"

/************************************************************************/
/*                         OGRMVTDBManager()                            */
/************************************************************************/
/**
 * @brief Constructs an OGRMVTDBManager object and prepare temp db path.
 * 
 * @param pszFilename The main dataset file name.
 * @param papszOptions Additional options containing configuration values.
 */
OGRMVTDBManager::OGRMVTDBManager(const char *pszFilename, char **papszOptions)
{
    m_osPath = GeneratePath(pszFilename, papszOptions);
}

/************************************************************************/
/*                        ~OGRMVTDBManager()                            */
/************************************************************************/
/**
 * @brief Destructor for OGRMVTDBManager class.
 * 
 * Finalizes all prepared SQLite statements and closes the database connection.
 * Ensures proper resource cleanup to prevent memory leaks and dangling database connection.
 */
OGRMVTDBManager::~OGRMVTDBManager()
{
    std::cout << "Destroying OGRMVTDBManager" << std::endl;
    CloseConnection();
}

/************************************************************************/
/*                          Initialize()                                */
/************************************************************************/
/**
 * @brief Initializes the temporary database (creates connection and table).
 * 
 * @param reuseExisting If true, it will reuse an existing database file.
 * @return OGRERR_NONE on success, or OGRERR_FAILURE on failure.
 */
OGRErr OGRMVTDBManager::Initialize(bool reuseExisting)
{
    std::cout << "OGRMVTDBManager::Initialize" << std::endl;

    if (OpenConnection(reuseExisting) != OGRERR_NONE)
        return OGRERR_FAILURE;

    return OGRERR_NONE;
}

/************************************************************************/
/*                             UnlinkFileIfNeeded                       */
/************************************************************************/
// @brief Method for unlinking (deleting) the temporary file
// based on the configuration and reuse settings.
void OGRMVTDBManager::UnlinkFileIfNeeded()
{
    std::cout << "OGRMVTDBManager::UnlinkFileIfNeeded" << std::endl;

    if (!m_osPath.empty() && !m_bReuseDB &&
        CPLTestBool(CPLGetConfigOption("OGR_MVT_REMOVE_TEMP_FILE", "YES")))
    {
        VSIUnlink(m_osPath);
    }
}

/************************************************************************/
/*                           CloseConnection()                          */
/************************************************************************/
/**
 * @brief Closes the database connection and finalizes all prepared statements.
 *
 * Ensures that all allocated SQLite statements are finalized and the database 
 * connection is properly closed. This prevents resource leaks and is called 
 * automatically in the destructor.
 */
void OGRMVTDBManager::CloseConnection()
{
    std::cout << "OGRMVTDBManager::CloseConnection" << std::endl;

    if (m_hInsertStmt)
        sqlite3_finalize(m_hInsertStmt);
    if (m_hTilesStmt)
        sqlite3_finalize(m_hTilesStmt);
    if (m_hLayersStmt)
        sqlite3_finalize(m_hLayersStmt);
    if (m_hFeaturesStmt)
        sqlite3_finalize(m_hFeaturesStmt);
    if (m_hFeatureLimitStmt)
        sqlite3_finalize(m_hFeatureLimitStmt);

    if (m_poDB)
    {
        sqlite3_close(m_poDB);
    }
}

/************************************************************************/
/*                              GeneratePath()                          */
/************************************************************************/
/**
 * @brief Generates a path for the temporary database.
 * 
 * @param pszFilename Path to the main dataset file.
 * @param papszOptions Additional options containing configuration values.
 * @return The generated path to the temporary database
 */
CPLString OGRMVTDBManager::GeneratePath(const char *pszFilename,
                                        char **papszOptions)
{
    std::cout << "OGRMVTDBManager::GeneratePath" << std::endl;

    CPLString osDBDefault = CPLString(pszFilename) + ".temp.db";
    if (STARTS_WITH(osDBDefault, "/vsizip/"))
    {
        osDBDefault = CPLString(pszFilename + strlen("/vsizip/")) + ".temp.db";
    }

    return CSLFetchNameValueDef(papszOptions, "TEMPORARY_DB",
                                osDBDefault.c_str());
}

/************************************************************************/
/*                          OpenConnection()                            */
/************************************************************************/
/**
 * @brief Opens or creates an SQLite database at the specified path.
 *
 * If the database file exists and `reuseExisting` is true, the method reuses it.
 * If `reuseExisting` is false, any existing file is deleted, and a new database is created.
 *
 * @param reuseExisting Whether to reuse an existing database file:
 *                       - true: Reuse the existing file if it exists.
 *                       - false: Delete the existing file and create a new one.
 *
 * @return OGRERR_NONE on success, or OGRERR_FAILURE on failure.
 */
OGRErr OGRMVTDBManager::OpenConnection(bool reuseExisting)
{
    std::cout << "OGRMVTDBManager::OpenConnection" << std::endl;

    bool dbExists = VSIStatL(m_osPath.c_str(), nullptr) == 0;

    // If not reusing and database exists, delete it
    if (!reuseExisting && dbExists)
    {
        VSIUnlink(m_osPath.c_str());
        dbExists = false;
    }

    int openFlags = dbExists ? SQLITE_OPEN_READWRITE | SQLITE_OPEN_NOMUTEX
                             : SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE |
                                   SQLITE_OPEN_NOMUTEX;

    if (sqlite3_open_v2(m_osPath.c_str(), &m_poDB, openFlags, nullptr) !=
            SQLITE_OK ||
        m_poDB == nullptr)
    {
        CPLError(CE_Failure, CPLE_FileIO, "Cannot open or create %s",
                 m_osPath.c_str());
        if (m_poDB)
        {
            sqlite3_close(m_poDB);
            m_poDB = nullptr;
        }
        return OGRERR_FAILURE;
    }

    CPLDebug("OGR",
             dbExists && reuseExisting ? "Reusing existing database: %s"
                                       : "Created new database: %s",
             m_osPath.c_str());

    return OGRERR_NONE;
}

/************************************************************************/
/*                         UpdateFeatureCount                       */
/************************************************************************/
// @brief Refresher for the number of features in the temporary database.
void OGRMVTDBManager::UpdateFeatureCount()
{
    std::cout << "OGRMVTDBManager::UpdateFeatureCount" << std::endl;

    const char *sqlQuery = "SELECT COUNT(*) FROM temp";
    m_nFeatureCount = SQLGetInteger64(m_poDB, sqlQuery, nullptr);
    ;
}

/************************************************************************/
/*                      CreateDataTable()                               */
/************************************************************************/
/**
 * Creates the temp table in the SQLite database used for storing MVT data.
 *
 * @return OGRERR_NONE on success, or OGRERR_FAILURE on failure.
 */
OGRErr OGRMVTDBManager::CreateDataTable()
{
    std::cout << "OGRMVTDBManager::CreateDataTable" << std::endl;

    const char *sql =
        "PRAGMA page_size = 4096;"  // 4096: default since SQLite 3.12
        "PRAGMA synchronous = OFF;"
        "PRAGMA journal_mode = OFF;"
        "PRAGMA temp_store = MEMORY;"
        "CREATE TABLE temp (z INTEGER, x INTEGER, y INTEGER, layer TEXT, "
        "idx INTEGER, feature BLOB, geomtype INTEGER, area_or_length DOUBLE);"
        "CREATE INDEX temp_index ON temp (z, x, y, layer, idx);";

    return SQLCommand(m_poDB, sql);
}

/************************************************************************/
/*                          PrepareInsertStmt()                             */
/************************************************************************/
/**
 * @brief Prepares an SQL statement to insert tiles into the temporary database.
 *
 * This method compiles an SQL query that inserts decoded tile features into
 * the `temp` table of the temporary SQLite database.
 *
 * @param hDB SQLite database handle for the temporary database.
 * @param stmt Pointer to the SQLite statement to be prepared.
 *
 * @return OGRERR_NONE on success, or OGRERR_FAILURE on failure.
 */
OGRErr OGRMVTDBManager::PrepareInsertFeatureStmt()
{
    std::cout << "OGRMVTDBManager::PrepareInsertStmt" << std::endl;

    const char *pszSQL = "INSERT INTO temp (z, x, y, layer, idx, feature, "
                         "geomtype, area_or_length) "
                         "VALUES (?, ?, ?, ?, ?, ?, ?, ?)";

    if (sqlite3_prepare_v2(m_poDB, pszSQL, -1, &m_hInsertStmt, nullptr) !=
        SQLITE_OK)
    {
        m_hInsertStmt = nullptr;
        return OGRERR_FAILURE;
    }
    return OGRERR_NONE;
}

/************************************************************************/
/*                          InsertFeature()                             */
/************************************************************************/
/**
 * @brief Inserts a single vector feature into the temp database.
 *
 * @param nZ        Zoom level of the tile.
 * @param nTileX    X coordinate of the tile.
 * @param nTileY    Y coordinate of the tile.
 * @param layerName Name of the target layer.
 * @param featureId Unique ID of the feature.
 * @param buffer    Serialized (protobuf) representation of the feature.
 * @param geomType  Geometry type (e.g. point, linestring, polygon).
 * @param areaOrLength Pre-computed area or length, if available.
 * 
 * @return OGRERR_NONE on success, OGRERR_FAILURE on error.
 */
OGRErr OGRMVTDBManager::InsertFeature(int nZ, int nTileX, int nTileY,
                                      const std::string &layerName,
                                      GIntBig featureId, CPLString buffer,
                                      int geomType, double areaOrLength)
{

    std::cout << "OGRMVTDBManager::InsertFeature" << std::endl;

    // Příprava příkazu INSERT, pokud ještě nebyl připraven
    if (m_hInsertStmt == nullptr)
    {
        PrepareInsertFeatureStmt();
        if (m_hInsertStmt == nullptr)
        {
            return OGRERR_FAILURE;  // Pokud příkaz není připraven, vrátí se chyba
        }
    }

    // Bindování hodnot a provedení SQL příkazu
    m_nFeatureCount++;

    sqlite3_bind_int(m_hInsertStmt, 1, nZ);
    sqlite3_bind_int(m_hInsertStmt, 2, nTileX);
    sqlite3_bind_int(m_hInsertStmt, 3, nTileY);
    sqlite3_bind_text(m_hInsertStmt, 4, layerName.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_int64(m_hInsertStmt, 5, featureId);
    sqlite3_bind_blob(m_hInsertStmt, 6, buffer.data(),
                      static_cast<int>(buffer.size()), SQLITE_STATIC);
    sqlite3_bind_int(m_hInsertStmt, 7, geomType);
    sqlite3_bind_double(m_hInsertStmt, 8, areaOrLength);

    int rc = sqlite3_step(m_hInsertStmt);
    if (rc != SQLITE_OK && rc != SQLITE_DONE)
    {
        sqlite3_reset(m_hInsertStmt);
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Failed to insert feature into temp DB: %s",
                 sqlite3_errmsg(m_poDB));
        return OGRERR_FAILURE;
    }

    sqlite3_reset(m_hInsertStmt);
    return OGRERR_NONE;
}

/************************************************************************/
/*                    FinalizeInsertFeatureStmt()                        */
/************************************************************************/
/**
 * @brief Finalizes the prepared statement for the feature insert query.
 */
void OGRMVTDBManager::FinalizeInsertFeatureStmt()
{
    if (m_hInsertStmt)
    {
        sqlite3_finalize(m_hInsertStmt);
        m_hInsertStmt = nullptr;
    }
}

/************************************************************************/
/*                        PrepareTilesStmt()                            */
/************************************************************************/
/**
 * @brief Prepares the SQL statement to fetch tiles in the dataset.
 * @return OGRERR_NONE if successful, OGRERR_FAILURE otherwise.
 */
OGRErr OGRMVTDBManager::PrepareTilesStmt()
{
    if (m_hTilesStmt == nullptr)
    {
        const char *sql = "SELECT DISTINCT z, x, y FROM temp ORDER BY z, x, y";
        if (sqlite3_prepare_v2(m_poDB, sql, -1, &m_hTilesStmt, nullptr) !=
            SQLITE_OK)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Failed to prepare m_hTilesStmt");
            return OGRERR_FAILURE;
        }
    }
    return OGRERR_NONE;
}

/************************************************************************/
/*                    FinalizeTilesStmt()                               */
/************************************************************************/
/**
 * @brief Finalizes the prepared statement for the tiles query.
 */
void OGRMVTDBManager::FinalizeTilesStmt()
{
    if (m_hTilesStmt)
    {
        sqlite3_finalize(m_hTilesStmt);
        m_hTilesStmt = nullptr;
    }
}

/************************************************************************/
/*                        PrepareLayersStmt                             */
/************************************************************************/
/**
 * @brief Prepares the SQL statement to fetch layers in a tile.
 * @return OGRERR_NONE if successful, OGRERR_FAILURE otherwise.
 */
OGRErr OGRMVTDBManager::PrepareLayersStmt()
{
    if (m_hLayersStmt == nullptr)
    {
        const char *sql = "SELECT DISTINCT layer FROM temp WHERE z = ? AND x = "
                          "? AND y = ? ORDER BY layer";
        if (sqlite3_prepare_v2(m_poDB, sql, -1, &m_hLayersStmt, nullptr) !=
            SQLITE_OK)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Failed to prepare m_hLayersStmt");
            FinalizeTilesStmt();
            return OGRERR_FAILURE;
        }
    }
    return OGRERR_NONE;
}

/************************************************************************/
/*                    FinalizeLayersStmt()                              */
/************************************************************************/
/**
 * @brief Finalizes the prepared statement for fetching layers in a tile.
 */
void OGRMVTDBManager::FinalizeLayersStmt()
{
    if (m_hLayersStmt)
    {
        sqlite3_finalize(m_hLayersStmt);
        m_hLayersStmt = nullptr;
    }
}

/************************************************************************/
/*                        PrepareFeaturesStmt                           */
/************************************************************************/
/**
 * @brief Prepares the SQL statement to fetch features in a given layer of a tile.
 * @return OGRERR_NONE if successful, OGRERR_FAILURE otherwise.
 */
OGRErr OGRMVTDBManager::PrepareFeaturesStmt()
{
    if (m_hFeaturesStmt == nullptr)
    {
        const char *sql = "SELECT feature FROM temp WHERE z = ? AND x = ? AND "
                          "y = ? AND layer = ? ORDER BY idx";
        if (sqlite3_prepare_v2(m_poDB, sql, -1, &m_hFeaturesStmt, nullptr) !=
            SQLITE_OK)
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Failed to prepare m_hFeaturesStmt");
            FinalizeTilesStmt();
            FinalizeLayersStmt();
            return OGRERR_FAILURE;
        }
    }
    return OGRERR_NONE;
}

/************************************************************************/
/*                    FinalizeFeaturesStmt()                            */
/************************************************************************/
/**
 * @brief Finalizes the prepared statement for fetching features in a given layer of a tile.
 */
void OGRMVTDBManager::FinalizeFeaturesStmt()
{
    if (m_hFeaturesStmt)
    {
        sqlite3_finalize(m_hFeaturesStmt);
        m_hFeaturesStmt = nullptr;
    }
}

/************************************************************************/
/*                      PrepareFeatureLimitStmt()                       */
/************************************************************************/
/**
 * @brief Prepares an SQLite statement to retrieve features from a temporary table.
 *
 * This function compiles an SQL query that selects layers and features from the 
 * `temp` table based on tile coordinates (z, x, y). The results are ordered by 
 * area or length in descending order and limited to a specified number of features.
 *
 * @return OGRERR_NONE on success, OGRERR_FAILURE on failure.
 */
OGRErr OGRMVTDBManager::PrepareFeatureLimitStmt()
{
    std::cout << "OGRMVTDBManager::PrepareFeatureLimitStmt" << std::endl;

    const char *pszSQL =
        "SELECT layer, feature FROM temp "
        "WHERE z = ? AND x = ? AND y = ? ORDER BY area_or_length DESC LIMIT ?";

    if (sqlite3_prepare_v2(m_poDB, pszSQL, -1, &m_hFeatureLimitStmt, nullptr) !=
        SQLITE_OK)
    {
        m_hFeatureLimitStmt = nullptr;
        return OGRERR_FAILURE;
    }
    return OGRERR_NONE;
}

/************************************************************************/
/*                    BindFeatureLimitStmtParams()                      */
/************************************************************************/
/**
 * @brief Binds parameters to the prepared SQL statement for limited feature retrieval.
 *
 * This method binds the tile coordinates (z, x, y) and a feature limit to the
 * previously prepared SQL statement (`m_hFeatureLimitStmt`). It also resets the 
 * statement to allow re-use with new parameters.
 *
 * @param nZ Zoom level.
 * @param nTileX Tile X coordinate.
 * @param nTileY Tile Y coordinate.
 * @param limit Maximum number of features to retrieve.
 *
 * @return OGRERR_NONE on success, OGRERR_FAILURE if binding fails or the statement is not prepared.
 */
OGRErr OGRMVTDBManager::BindFeatureLimitStmtParams(int nZ, int nTileX,
                                                   int nTileY, unsigned limit)
{
    if (!m_hFeatureLimitStmt)
        return OGRERR_FAILURE;

    sqlite3_reset(m_hFeatureLimitStmt);  // umožní opakované použití

    if (sqlite3_bind_int(m_hFeatureLimitStmt, 1, nZ) != SQLITE_OK ||
        sqlite3_bind_int(m_hFeatureLimitStmt, 2, nTileX) != SQLITE_OK ||
        sqlite3_bind_int(m_hFeatureLimitStmt, 3, nTileY) != SQLITE_OK ||
        sqlite3_bind_int(m_hFeatureLimitStmt, 4, static_cast<int>(limit)) !=
            SQLITE_OK)
    {
        return OGRERR_FAILURE;
    }
    return OGRERR_NONE;
}

/************************************************************************/
/*                    FinalizeFeatureLimitStmt()                        */
/************************************************************************/
/**
 * @brief Finalizes the prepared statement for the feature limit query.
 */
void OGRMVTDBManager::FinalizeFeatureLimitStmt()
{
    if (m_hFeatureLimitStmt)
    {
        sqlite3_finalize(m_hFeatureLimitStmt);
        m_hFeatureLimitStmt = nullptr;
    }
}

/************************************************************************/
/*                        PrepareOutputStmts                            */
/************************************************************************/
/**
 * @brief Prepares all required SQL statements for output queries.
 * @return OGRERR_NONE if all statements are prepared successfully, OGRERR_FAILURE otherwise.
 */
OGRErr OGRMVTDBManager::PrepareOutputStmts()
{
    if (PrepareTilesStmt() != OGRERR_NONE)
        return OGRERR_FAILURE;
    if (PrepareLayersStmt() != OGRERR_NONE)
        return OGRERR_FAILURE;
    if (PrepareFeaturesStmt() != OGRERR_NONE)
        return OGRERR_FAILURE;
    return OGRERR_NONE;
}

/************************************************************************/
/*                        FinalizeOutputStmts()                        */
/************************************************************************/
/**
 * @brief Finalizes output SQL statements used for retrieving tile-related data.
 */
void OGRMVTDBManager::FinalizeOutputStmts()
{
    FinalizeTilesStmt();
    FinalizeLayersStmt();
    FinalizeFeaturesStmt();
}