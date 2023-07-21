#include "hdf_base_depth_reader.h"

//------------------------------------------------------------------------------
bool HdfBaseDepthReader::open(const char *filename) {
    if (isFileOpened) {
        delete[] hdfFile;
        isFileOpened = false;
    }

    try {
        H5::Exception::dontPrint();
        hdfFile = new H5::H5File(filename, H5F_ACC_RDONLY);
    } catch (H5::FileIException err) {
        std::cerr << "The specified HDF matrix (" << filename << ") can't open."
                  << ENDL;
        return false;
    }

    isFileOpened = true;
    return true;
}
//------------------------------------------------------------------------------
bool HdfBaseDepthReader::set_target_chromosome(const char *chr_str) {
    if (chr_str != currentChrName) {
        if (isGroupAllocated) {
            delete hdfGroup;
            isGroupAllocated = false;
        }

        try {
            H5::Exception::dontPrint();
            hdfGroup = new H5::Group(hdfFile->openGroup(chr_str));
            isGroupAllocated = true;
        } catch (H5::GroupIException err) {
            std::cerr << ERROR_STRING << "a replicon '" << chr_str
                      << "' can't open." << ENDL;
            return false;
        }
        currentChrName = chr_str;
    }
    return true;
}
//------------------------------------------------------------------------------
bool HdfBaseDepthReader::set_target_dataset(
        const char *dataName, const H5::DataType dataType) {
    if (isDataSetAllocated) {
        delete hdfDataSet;
        isDataSetAllocated = false;
    }
    if (isDataSpaceAllocated) {
        delete hdfDataSpace;
        isDataSpaceAllocated = false;
    }

    try {
        hdfDataSet = new H5::DataSet(hdfGroup->openDataSet(dataName));
        isDataSetAllocated = true;
        hdfDataSpace = new H5::DataSpace(hdfDataSet->getSpace());
        isDataSpaceAllocated = true;
    } catch (H5::GroupIException err) {
        std::cerr << ERROR_STRING << "can't find the specified dataset '"
                  << dataName << "' in the group." << ENDL;
        return false;
    } catch (H5::DataSetIException err) {
        std::cerr << ERROR_STRING << "can't find a specified dataset '"
                  << dataName << "'." << ENDL;
        return false;
    } catch (H5::DataSpaceIException err) {
        std::cerr << ERROR_STRING << "can't create a dataspace." << ENDL;
        return false;
    }

    currentDataName = dataName;
    currentDataType = dataType;
    return true;
}
//------------------------------------------------------------------------------
bool HdfBaseDepthReader::get_matrix(
        const int *start, const int *count, IntType *buffer) {
    if (!isFileOpened)
        return false;

    hsize_t f_offset[1], m_offset[1], h_count[1];
    f_offset[0] = *start;
    m_offset[0] = 0;
    h_count[0] = *count;

    try {
        H5::Exception::dontPrint();
        H5::DataSpace *memspace = new H5::DataSpace(1, h_count);
        memspace->selectHyperslab(H5S_SELECT_SET, h_count, m_offset);
        hdfDataSpace->selectHyperslab(H5S_SELECT_SET, h_count, f_offset);
        hdfDataSet->read(
                (void *)buffer, currentDataType, *memspace, *hdfDataSpace);
        delete memspace;
    } catch (H5::DataSetIException err) {
        std::cerr << ERROR_STRING << "can't read data from a dataspace."
                  << " start=" << *start << ", count=" << *count << ENDL;
        return false;
    }
    return true;
}
//------------------------------------------------------------------------------
bool HdfBaseDepthReader::get_group_names(hi::StringArray &groups) {
    H5::Group *group = new H5::Group(hdfFile->openGroup("/"));
    hi::StringArray *ptr = &groups;
    if (0
            == H5Literate(group->getId(), H5_INDEX_NAME, H5_ITER_INC, NULL,
                    add_group, (void *)ptr))
        return true;
    return false;
}
//------------------------------------------------------------------------------
bool HdfBaseDepthReader::get_unique_read_count(
        const char *chr, int *read_count) {
    if (set_target_chromosome(chr)
            && set_target_dataset("UniqueReadCount", H5::PredType::STD_I32LE)) {
        hdfDataSet->read(read_count, H5::PredType::STD_I32LE, *hdfDataSpace);
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
int HdfBaseDepthReader::get_num_elements(void) {
    if (!isDataSpaceAllocated)
        return -1;

    hsize_t dims[2];
    int ndims = hdfDataSpace->getSimpleExtentDims(dims, NULL);
    if (0 > ndims)
        return -1;

    return dims[0];
}
//------------------------------------------------------------------------------
void HdfBaseDepthReader::close(void) {
    if (isFileOpened) {
        delete hdfFile;
        isFileOpened = false;
    }
    if (isDataSpaceAllocated) {
        delete hdfDataSpace;
        isDataSpaceAllocated = false;
    }
    if (isDataSetAllocated) {
        delete hdfDataSet;
        isDataSetAllocated = false;
    }
    if (isGroupAllocated) {
        delete hdfGroup;
        isGroupAllocated = false;
    }
}
//------------------------------------------------------------------------------
HdfBaseDepthReader::HdfBaseDepthReader() {
    isFileOpened = false;
    isGroupAllocated = false;
    isDataSetAllocated = false;
    isDataSpaceAllocated = false;
}
//------------------------------------------------------------------------------
HdfBaseDepthReader::~HdfBaseDepthReader() {
    this->close();
}
//------------------------------------------------------------------------------
herr_t add_group(hid_t loc_id, const char *namestr, const H5L_info_t *linfo,
        void *retval) {
    hi::StringArray *names = (hi::StringArray *)retval;
    names->push_back(namestr);
    return 0;
}
//------------------------------------------------------------------------------
