/******************************************************************************
 *
 * Project:  MVT Translator
 * Purpose:  Mapbox Vector Tile decoder and encoder
 * Author:   Even Rouault, Even Rouault <even dot rouault at spatialys dot com>
 *           Linda Karlovska <linda dot karlovska at seznam dot cz> (refactoring)
 *
 ******************************************************************************
 * Copyright (c) 2018, Even Rouault <even dot rouault at spatialys dot com>
 *               2025, Linda Karlovska <linda dot karlovska at seznam dot cz>
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#ifndef OGRMVTDB_H
#define OGRMVTDB_H

#include <mutex>
#include <vector>
#include <tuple>
#include <iostream>
#include <sqlite3.h>
#include "../sqlite/ogrsqliteutility.h"

/************************************************************************/
/*                         OGRMVTDBManager                              */
/************************************************************************/
/** 
 * \brief Helper class managing temporary SQLite database used for
 *        intermediate feature storage during MVT writing process.
 *
 * This class provides methods for creating the schema, preparing SQL
 * statements, and inserting or querying feature-level data during the 
 * vector tile generation process. 
 */
class OGRMVTDBManager
{
  public:
    explicit OGRMVTDBManager(const char *pszFilename, char **papszOptions);
    ~OGRMVTDBManager();

    OGRErr Initialize(bool reuseExisting);
    void UnlinkFileIfNeeded();
    void CloseConnection();

    void UpdateFeatureCount();

    OGRErr CreateDataTable();

    CPLString GetPath() const
    {
        return m_osPath;
    }

    sqlite3 *GetDB() const
    {
        return m_poDB;
    }

    bool GetReuseDB() const
    {
        return m_bReuseDB;
    }

    GIntBig GetFeatureCount() const
    {
        return m_nFeatureCount;
    }

    // --- Prepare statements ---
    OGRErr PrepareInsertFeatureStmt();
    OGRErr PrepareTilesStmt();
    OGRErr PrepareLayersStmt();
    OGRErr PrepareFeaturesStmt();
    OGRErr PrepareFeatureLimitStmt();
    OGRErr BindFeatureLimitStmtParams(int nZ, int nTileX, int nTileY,
                                      unsigned limit);
    OGRErr PrepareOutputStmts();

    // --- Getters for prepared statements ---
    sqlite3_stmt *GetInsertFeatureStmt() const
    {
        return m_hInsertStmt;
    }

    sqlite3_stmt *GetTilesStmt() const
    {
        return m_hTilesStmt;
    }

    sqlite3_stmt *GetLayersStmt() const
    {
        return m_hLayersStmt;
    }

    sqlite3_stmt *GetFeaturesStmt() const
    {
        return m_hFeaturesStmt;
    }

    sqlite3_stmt *GetFeatureLimitStmt() const
    {
        return m_hFeatureLimitStmt;
    }

    // --- Finalize statements ---
    void FinalizeInsertFeatureStmt();
    void FinalizeTilesStmt();
    void FinalizeLayersStmt();
    void FinalizeFeaturesStmt();
    void FinalizeFeatureLimitStmt();
    void FinalizeOutputStmts();

    // --- Feature writing ---
    OGRErr InsertFeature(int nZ, int nTileX, int nTileY,
                         const std::string &layerName, GIntBig featureId,
                         CPLString buffer, int geomType, double areaOrLength);

  private:
    CPLString m_osPath;
    sqlite3 *m_poDB = nullptr;
    mutable GIntBig m_nFeatureCount = 0;
    bool m_bReuseDB = false;

    // --- Prepared SQLite statements ---
    sqlite3_stmt *m_hInsertStmt = nullptr;
    sqlite3_stmt *m_hTilesStmt = nullptr;
    sqlite3_stmt *m_hLayersStmt = nullptr;
    sqlite3_stmt *m_hFeaturesStmt = nullptr;
    sqlite3_stmt *m_hFeatureLimitStmt = nullptr;

    CPLString GeneratePath(const char *pszFilename, char **papszOptions);
    OGRErr OpenConnection(bool reuseExisting);
};

#endif  // OGRMVTDB_H
