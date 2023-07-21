#ifndef HDF_BASE_DEPTH_READER_H
#define HDF_BASE_DEPTH_READER_H

#include "histd.h"

#include <H5Cpp.h>
#include <iostream>
#include <string>

//------------------------------------------------------------------------------
typedef int32_t IntType;
//------------------------------------------------------------------------------
class HdfBaseDepthReader {
public:
    bool open(const char *filename);
    bool set_target_chromosome(const char *chr_str);
    bool set_target_dataset(const char *dataName, const H5::DataType dataType);
    bool get_matrix(const int *start, const int *count, IntType *buffer);
    bool get_group_names(hi::StringArray &groups);
    bool get_unique_read_count(const char *chr, int *read_count);
    int get_num_elements(void);
    void close(void);
    HdfBaseDepthReader();
    ~HdfBaseDepthReader();

protected:
    H5::H5File *hdfFile;
    H5::Group *hdfGroup;
    H5::DataSet *hdfDataSet;
    H5::DataSpace *hdfDataSpace;

    std::string currentChrName, currentDataName;
    H5::DataType currentDataType;
    bool isFileOpened, isGroupAllocated, isDataSetAllocated,
            isDataSpaceAllocated;
};
//------------------------------------------------------------------------------
herr_t add_group(hid_t loc_id, const char *namestr, const H5L_info_t *linfo,
        void *retval);
//------------------------------------------------------------------------------

#endif
