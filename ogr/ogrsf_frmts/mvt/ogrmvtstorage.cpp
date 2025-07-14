/******************************************************************************
 *
 * Project:  MVT Translator
 * Purpose:  MVT storage management implementation, including factory.
 * Author:   Even Rouault, Even Rouault <even dot rouault at spatialys dot com>
 *           Linda Karlovska <linda dot karlovska at seznam dot cz>
 *
 ******************************************************************************
 * Copyright (c) 2018, Even Rouault <even dot rouault at spatialys dot com>
 * Copyright (c) 2025, Linda Karlovska <linda dot karlovska at seznam dot cz>
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#include "ogrmvtstorage.h"
#include "ogrsf_frmts.h"
#include "cpl_vsi.h"
#include "cpl_string.h"
#include "cpl_error.h"

#include <sqlite3.h>
#include "../sqlite/ogrsqliteutility.h"

#include <cmath>

/************************************************************************/
/*                            CreateStorageManager()                    */
/************************************************************************/
std::unique_ptr<OGRMVTStorageManager>
CreateStorageManager(const char *pszFilename, char **papszOptions)
{
    const char *pszFormat = CSLFetchNameValue(papszOptions, "FORMAT");
    const bool bMBTILESExt =
        EQUAL(CPLGetExtensionSafe(pszFilename).c_str(), "mbtiles");
    if (pszFormat == nullptr && bMBTILESExt)
        pszFormat = "MBTILES";
    const bool bMBTILES = pszFormat != nullptr && EQUAL(pszFormat, "MBTILES");

    const CPLString osFilename(pszFilename);
    const CPLString osExtension =
        CSLFetchNameValueDef(papszOptions, "TILE_EXTENSION", "");

    if (bMBTILES)
        return std::make_unique<MBTilesStorageManager>(pszFilename,
                                                       osExtension);
    else
        return std::make_unique<DirectoryStorageManager>(pszFilename,
                                                         osExtension);
}

/************************************************************************/
/*                         DirectoryStorageManager()                    */
/************************************************************************/
/************************************************************************/
/*                             Initialize()                             */
/************************************************************************/
bool DirectoryStorageManager::Initialize(CPL_UNUSED const char *pszVFSName)
{
    VSIStatBufL sStat;
    const char *osFilename = GetFilename().c_str();

    if (VSIStatL(osFilename, &sStat) == 0)
    {
        CPLError(CE_Failure, CPLE_FileIO, "%s already exists", osFilename);
        return false;
    }

    if (VSIMkdir(osFilename, 0755) != 0)
    {
        CPLError(CE_Failure, CPLE_FileIO, "Cannot create directory %s",
                 osFilename);
        return false;
    }

    return true;
}

/************************************************************************/
/*                              WriteTile()                             */
/************************************************************************/
bool DirectoryStorageManager::WriteTile(const std::string &oTileBuffer, int nZ,
                                        int nX, int nY)
{
    bool bRet = true;

    CPLString osZDirname(CPLFormFilenameSafe(GetFilename().c_str(),
                                             CPLSPrintf("%d", nZ), nullptr));
    CPLString osXDirname(
        CPLFormFilenameSafe(osZDirname, CPLSPrintf("%d", nX), nullptr));

    // Create directories if they don't exist
    if (nZ != m_nLastZ)
    {
        VSIMkdir(osZDirname, 0755);  // Create directory for zoom level
        m_nLastZ = nZ;
        m_nLastX = -1;
    }
    if (nX != m_nLastX)
    {
        VSIMkdir(osXDirname, 0755);  // Create directory for X coordinate
        m_nLastX = nX;
    }

    // Generate the tile filename
    CPLString osTileFilename(CPLFormFilenameSafe(
        osXDirname, CPLSPrintf("%d", nY), GetExtension().c_str()));

    // Open the file for writing
    VSILFILE *fpOut = VSIFOpenL(osTileFilename, "wb");
    if (fpOut)
    {
        // Write the tile buffer to the file
        const size_t nRet =
            VSIFWriteL(oTileBuffer.data(), 1, oTileBuffer.size(), fpOut);
        bRet = (nRet == oTileBuffer.size());  // Check if write was successful
        VSIFCloseL(fpOut);                    // Close the file after writing
    }
    else
    {
        bRet = false;  // Failed to open file
    }
    return bRet;
}

/************************************************************************/
/*                           WriteMetadataItem()                        */
/************************************************************************/
bool DirectoryStorageManager::WriteMetadataItem(const char *pszKey,
                                                const char *pszValue)
{
    return WriteMetadataItemInternal(pszKey, pszValue);
}

/************************************************************************/
/*                           WriteMetadataItem()                        */
/************************************************************************/
bool DirectoryStorageManager::WriteMetadataItem(const char *pszKey, int nValue)
{
    return WriteMetadataItemInternal(pszKey, nValue);
}

/************************************************************************/
/*                           WriteMetadataItem()                        */
/************************************************************************/
bool DirectoryStorageManager::WriteMetadataItem(const char *pszKey,
                                                double dfValue)
{
    return WriteMetadataItemInternal(pszKey, dfValue);
}

/************************************************************************/
/*                           WriteMetadataItemT()                       */
/************************************************************************/
template <typename T>
bool DirectoryStorageManager::WriteMetadataItemInternal(const char *pszKey,
                                                        T value)
{
    m_oRoot.Add(pszKey, value);
    return true;
}

/************************************************************************/
/*                           SaveMetadata()                             */
/************************************************************************/
bool DirectoryStorageManager::SaveMetadata()
{
    return m_oDoc.Save(CPLFormFilenameSafe(GetFilename(), "metadata.json",
                                           nullptr)) == OGRERR_NONE;
}

/************************************************************************/
/*                              Close()                                 */
/************************************************************************/
void DirectoryStorageManager::Close()
{
    // No-op: nothing to close in directory-based storage
}

/************************************************************************/
/*                            MBTilesStorageManager()                   */
/************************************************************************/

/************************************************************************/
/*                           ~MBTilesStorageManager()                   */
/************************************************************************/
MBTilesStorageManager::~MBTilesStorageManager()
{
    Close();
}

/************************************************************************/
/*                           Initialize()                    		    */
/************************************************************************/
bool MBTilesStorageManager::Initialize(const char *pszVFSName)
{
    const CPLString osFilenameStr = GetFilename();
    const char *osFilename = osFilenameStr.c_str();

    const bool bMBTILESExt =
        EQUAL(CPLGetExtensionSafe(osFilename).c_str(), "mbtiles");

    if (!bMBTILESExt)
    {
        CPLError(CE_Failure, CPLE_FileIO, "%s should have mbtiles extension",
                 osFilename);
        return false;
    }

    VSIUnlink(osFilename);

    if (sqlite3_open_v2(osFilename, &m_hDB,
                        SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE |
                            SQLITE_OPEN_NOMUTEX,
                        pszVFSName) != SQLITE_OK ||
        m_hDB == nullptr)
    {
        CPLError(CE_Failure, CPLE_FileIO, "Cannot create %s", osFilename);
        return false;
    }

    if (SQLCommand(
            m_hDB,
            "PRAGMA page_size = 4096;"  // 4096: default since sqlite 3.12
            "PRAGMA synchronous = OFF;"
            "PRAGMA journal_mode = OFF;"
            "PRAGMA temp_store = MEMORY;"
            "CREATE TABLE metadata (name text, value text);"
            "CREATE TABLE tiles (zoom_level integer, tile_column integer, "
            "tile_row integer, tile_data blob, "
            "UNIQUE (zoom_level, tile_column, tile_row))") != OGRERR_NONE)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Failed to initialize MBTiles DB");
        return false;
    }
    return true;
}

/************************************************************************/
/*                               WriteTile()                            */
/************************************************************************/
bool MBTilesStorageManager::WriteTile(const std::string &oTileBuffer, int nZ,
                                      int nX, int nY)
{
    sqlite3_stmt *hStmt = nullptr;
    const char *pszSQL = "INSERT INTO tiles(zoom_level, tile_column, tile_row, "
                         "tile_data) VALUES (?,?,?,?)";

    if (sqlite3_prepare_v2(m_hDB, pszSQL, -1, &hStmt, nullptr) != SQLITE_OK)
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Failed to prepare MBTiles insert");
        return false;
    }

    sqlite3_bind_int(hStmt, 1, nZ);
    sqlite3_bind_int(hStmt, 2, nX);
    sqlite3_bind_int(hStmt, 3, (1 << nZ) - 1 - nY);  // Flip Y axis
    sqlite3_bind_blob(hStmt, 4, oTileBuffer.data(),
                      static_cast<int>(oTileBuffer.size()), SQLITE_STATIC);

    const int rc = sqlite3_step(hStmt);
    bool bRet = (rc == SQLITE_OK || rc == SQLITE_DONE);

    sqlite3_finalize(hStmt);
    return bRet;
}

/************************************************************************/
/*                           WriteMetadataItem()                        */
/************************************************************************/
bool MBTilesStorageManager::WriteMetadataItem(const char *pszKey,
                                              const char *pszValue)
{
    return WriteMetadataItemInternal(pszKey, pszValue, "%q");
}

/************************************************************************/
/*                           WriteMetadataItem()                        */
/************************************************************************/
bool MBTilesStorageManager::WriteMetadataItem(const char *pszKey, int nValue)
{
    return WriteMetadataItemInternal(pszKey, nValue, "%d");
}

/************************************************************************/
/*                           WriteMetadataItem()                        */
/************************************************************************/
bool MBTilesStorageManager::WriteMetadataItem(const char *pszKey,
                                              double dfValue)
{
    return WriteMetadataItemInternal(pszKey, dfValue, "%.17g");
}

/************************************************************************/
/*                           WriteMetadataItemT                         */
/************************************************************************/
template <typename T>
bool MBTilesStorageManager::WriteMetadataItemInternal(
    const char *pszKey, T value, const char *pszValueFormat)
{
    char *pszSQL =
        sqlite3_mprintf("INSERT INTO metadata(name, value) VALUES(%Q, '%s')",
                        pszKey, CPLSPrintf(pszValueFormat, value));

    OGRErr eErr = SQLCommand(m_hDB, pszSQL);
    sqlite3_free(pszSQL);

    return eErr == OGRERR_NONE;
}

/************************************************************************/
/*                           SaveMetadata                               */
/************************************************************************/
bool MBTilesStorageManager::SaveMetadata()
{
    // No-op: nothing to save, metadata are already stored in the MBTiles DB
    return true;
}

/************************************************************************/
/*                              Close()                                 */
/************************************************************************/
void MBTilesStorageManager::Close()
{
    if (m_hDB)
    {
        sqlite3_close(m_hDB);
        m_hDB = nullptr;
    }
}