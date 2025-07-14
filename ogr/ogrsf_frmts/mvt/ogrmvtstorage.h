/******************************************************************************
 *
 * Project:  MVT Translator
 * Purpose:  MVT storage management base class and factory function.
 * Author:   Even Rouault, Even Rouault <even dot rouault at spatialys dot com>
 *           Linda Karlovska <linda dot karlovska at seznam dot cz>
 *
 ******************************************************************************
 * Copyright (c) 2018, Even Rouault <even dot rouault at spatialys dot com>
 * Copyright (c) 2025, Linda Karlovska <linda dot karlovska at seznam dot cz>
 * SPDX-License-Identifier: MIT
 ****************************************************************************/

#ifndef OGRMVTSTORAGE_H
#define OGRMVTSTORAGE_H

#include "cpl_json.h"

struct sqlite3;

/************************************************************************/
/*                         OGRMVTStorageManager                         */
/************************************************************************/
/**
 * \brief Abstract base class for storage of MVT tiles and metadata.
 *
 * Provides a common interface for writing tiles and associated metadata
 * into various storage formats (e.g., directory structure, MBTiles).
 */
class OGRMVTStorageManager
{
  public:
    explicit OGRMVTStorageManager(const CPLString &pszFilename,
                                  const CPLString &pszExtension)
        : m_osFilename(pszFilename),
          m_osExtension(pszExtension.empty() ? "pbf" : pszExtension)
    {
    }

    virtual ~OGRMVTStorageManager() = default;

    virtual bool Initialize(const char *pszVFSName = nullptr) = 0;

    virtual bool WriteTile(const std::string &oTileBuffer, int nZ, int nX,
                           int nY) = 0;

    virtual bool WriteMetadataItem(const char *pszKey,
                                   const char *pszValue) = 0;
    virtual bool WriteMetadataItem(const char *pszKey, int nValue) = 0;
    virtual bool WriteMetadataItem(const char *pszKey, double dfValue) = 0;

    virtual bool SaveMetadata() = 0;

    virtual bool IsMBTiles() const = 0;

    virtual void Close() = 0;

    const CPLString &GetFilename() const
    {
        return m_osFilename;
    }

    const CPLString &GetExtension() const
    {
        return m_osExtension;
    }

  protected:
    CPLString m_osFilename;
    CPLString m_osExtension{"pbf"};

  private:
    // Prevent copying
    OGRMVTStorageManager(const OGRMVTStorageManager &) = delete;
    OGRMVTStorageManager &operator=(const OGRMVTStorageManager &) = delete;
};

/************************************************************************/
/*                          CreateStorageManager()                      */
/************************************************************************/
/**
 * \brief Factory function to create an appropriate storage manager instance.
 */
std::unique_ptr<OGRMVTStorageManager>
CreateStorageManager(const char *pszFilename, char **papszOptions);

/************************************************************************/
/*                           DirectoryStorageManager()                  */
/************************************************************************/
/**
 * \brief Stores tiles in a directory structure and writes metadata as JSON.
 */
class DirectoryStorageManager final : public OGRMVTStorageManager
{
  public:
    DirectoryStorageManager(const CPLString &pszFilename,
                            const CPLString &pszExtension)
        : OGRMVTStorageManager(pszFilename, pszExtension)
    {
        m_oRoot = m_oDoc.GetRoot();
    }

    bool Initialize(CPL_UNUSED const char *pszVFSName) override;

    bool WriteTile(const std::string &oTileBuffer, int nZ, int nX,
                   int nY) override;

    bool WriteMetadataItem(const char *pszKey, const char *pszValue) override;
    bool WriteMetadataItem(const char *pszKey, int nValue) override;
    bool WriteMetadataItem(const char *pszKey, double dfValue) override;

    bool SaveMetadata() override;

    bool IsMBTiles() const override
    {
        return false;
    }

    void Close() override;

  private:
    CPLJSONDocument m_oDoc;
    CPLJSONObject m_oRoot;

    int m_nLastZ = -1;
    int m_nLastX = -1;

    template <typename T>
    bool WriteMetadataItemInternal(const char *pszKey, T value);
};

/************************************************************************/
/*                        MBTilesStorageManager                         */
/************************************************************************/
/**
 * \brief Stores tiles and metadata in an MBTiles SQLite database.
 */
class MBTilesStorageManager final : public OGRMVTStorageManager
{
  public:
    MBTilesStorageManager(const CPLString &pszFilename,
                          const CPLString &pszExtension)
        : OGRMVTStorageManager(pszFilename, pszExtension)
    {
    }

    ~MBTilesStorageManager() override;

    bool Initialize(const char *pszVFSName) override;

    bool WriteTile(const std::string &oTileBuffer, int nZ, int nX,
                   int nY) override;

    bool WriteMetadataItem(const char *pszKey, const char *pszValue) override;
    bool WriteMetadataItem(const char *pszKey, int nValue) override;
    bool WriteMetadataItem(const char *pszKey, double dfValue) override;

    bool SaveMetadata() override;

    bool IsMBTiles() const override
    {
        return true;
    }

    void Close() override;

  private:
    sqlite3 *m_hDB = nullptr;

    template <typename T>
    bool WriteMetadataItemInternal(const char *pszKey, T value,
                                   const char *pszValueFormat);
};

#endif  //OGRMVTSTORAGE_H