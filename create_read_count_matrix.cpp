#include "histd.h"

#include <H5Cpp.h>
#include <algorithm>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <malloc.h>
#include <sstream>
#include <string>
#include <unistd.h>

#define MAX_NAME_SIZE 255
#define NUM_THREADS 8
#define CHUNK_SIZE 65536

typedef int32_t IntMatrixType;
//------------------------------------------------------------------------------
struct ThreadCountParam {
    std::string input_fn, ref_name;
    int ref_start, ref_end, ref_length, unique_read_count;
    IntMatrixType *matrix, *clipend_matrix;
};
//------------------------------------------------------------------------------
void *thread_make_matrix(void *arg) {
    ThreadCountParam *p = (ThreadCountParam *)arg;

    // open the input bam & index files
    BamTools::BamReader bam_reader;
    if (!bam_reader.Open(p->input_fn)) {
        std::cerr << ERROR_STRING << "bam_reader.Open() failed at line "
                  << __LINE__ << ". input_fn=" << p->input_fn << ENDL;
        return NULL;
    }
    if (!bam_reader.LocateIndex(BamTools::BamIndex::STANDARD))
        bam_reader.CreateIndex();

    // set target region
    int refid = bam_reader.GetReferenceID(p->ref_name);
    bam_reader.SetRegion(refid, p->ref_start, refid, p->ref_end);

    // iterate through all alignments
    int alignment_start, alignment_end, endpos;
    BamTools::BamAlignment alignment;

    while (bam_reader.GetNextAlignmentCore(alignment)) {
        if (!alignment.IsMapped())
            continue;
        endpos = alignment.GetEndPosition();
        alignment_start = std::max(0, std::min(alignment.Position, endpos) - 1);
        alignment_end
                = std::min(p->ref_length, std::max(alignment.Position, endpos));
        for (int i = alignment_start; i < alignment_end; ++i)
            ++(p->matrix[i]);

        // soft-clipped read ends
        if ('S' == alignment.CigarData.at(0).Type
                && p->ref_length > alignment.Position && 0 <= endpos)
            p->clipend_matrix[alignment.Position]
                    += alignment.CigarData.at(0).Length;
        if ('S' == alignment.CigarData.at(alignment.CigarData.size() - 1).Type
                && p->ref_length > endpos && 0 <= endpos)
            p->clipend_matrix[endpos]
                    += alignment.CigarData.at(alignment.CigarData.size() - 1)
                               .Length;

        // count unique reads
        if (alignment.IsPrimaryAlignment())
            ++(p->unique_read_count);
    }

    bam_reader.Close();
    return NULL;
}
//------------------------------------------------------------------------------
BamTools::RefVector get_refvector(const std::string &input_fn) {
    // open the input bam & index files
    BamTools::BamReader bam_reader;
    if (!bam_reader.Open(input_fn)) {
        std::cerr << ERROR_STRING << "bam_reader.Open() failed at line "
                  << __LINE__ << ". input_fn=" << input_fn << ENDL;
        BamTools::RefVector dummy;
        return dummy;
    }
    if (!bam_reader.LocateIndex(BamTools::BamIndex::STANDARD))
        bam_reader.CreateIndex();

    BamTools::RefVector refvector = bam_reader.GetReferenceData();
    bam_reader.Close();

    return refvector;
}
//------------------------------------------------------------------------------
void write_hdf(H5::H5File *file, const char *ref_name,
        const IntMatrixType *matrix, const IntMatrixType *clipend_matrix,
        const int *szMatrix, const int *unique_read_count) {
    int rank = 1;
    hsize_t dims[2], cdims[2];
    H5::Group *group;
    H5::DataSet *dataset;

    // group for depth matrix
    std::stringstream fstr;
    fstr << "/" << ref_name;
    try {
        H5::Exception::dontPrint();
        group = new H5::Group(file->createGroup(fstr.str().c_str()));
    } catch (H5::GroupIException err) {
        std::cerr << ERROR_STRING << "group '" << fstr.str() << "' can't open."
                  << ENDL;
        return;
    }

    // write name
    fstr << "/FullName";
    dims[0] = std::strlen(ref_name);
    H5::DataSpace dataspace(rank, dims);
    try {
        H5::Exception::dontPrint();
        dataset = new H5::DataSet(file->createDataSet(
                fstr.str().c_str(), H5::PredType::C_S1, dataspace));
    } catch (H5::DataSetIException err) {
        std::cerr << ERROR_STRING << "dataset '" << fstr.str()
                  << "' can't open. func_name='" << err.getFuncName()
                  << "', msg='" << err.getDetailMsg() << "'." << ENDL;
        return;
    }
    dataset->write(ref_name, H5::PredType::C_S1, dataspace);
    delete dataset;

    // write length
    fstr.str("/");
    fstr << ref_name << "/Length";
    dims[0] = 1;
    H5::DataSpace dataspace2(rank, dims);
    try {
        H5::Exception::dontPrint();
        dataset = new H5::DataSet(file->createDataSet(
                fstr.str().c_str(), H5::PredType::STD_I32LE, dataspace2));
    } catch (H5::DataSetIException err) {
        std::cerr << ERROR_STRING << "dataset '" << fstr.str()
                  << "' can't open. func_name='" << err.getFuncName()
                  << "', msg='" << err.getDetailMsg() << "'." << ENDL;
        return;
    }
    dataset->write(szMatrix, H5::PredType::STD_I32LE, dataspace2);
    delete dataset;

    // unique read count -- Added in the matrix Version 0.2
    fstr.str("/");
    fstr << ref_name << "/UniqueReadCount";
    try {
        H5::Exception::dontPrint();
        dataset = new H5::DataSet(file->createDataSet(
                fstr.str().c_str(), H5::PredType::STD_I32LE, dataspace2));
    } catch (H5::DataSetIException err) {
        std::cerr << ERROR_STRING << "dataset '" << fstr.str()
                  << "' can't open. func_name='" << err.getFuncName()
                  << "', msg='" << err.getDetailMsg() << "'." << ENDL;
        return;
    }
    dataset->write(unique_read_count, H5::PredType::STD_I32LE, dataspace2);
    delete dataset;

    // base depth
    fstr.str("/");
    fstr << ref_name << "/BaseDepth";
    dims[0] = *szMatrix;
    cdims[0] = std::min(CHUNK_SIZE, *szMatrix);
    H5::DataSpace dataspace3(rank, dims);

    H5::DSetCreatPropList ds_creatplist;
    ds_creatplist.setChunk(rank, cdims);
    ds_creatplist.setDeflate(5);

    try {
        H5::Exception::dontPrint();
        dataset = new H5::DataSet(file->createDataSet(fstr.str().c_str(),
                H5::PredType::STD_I32LE, dataspace3, ds_creatplist));
    } catch (H5::DataSetIException err) {
        std::cerr << ERROR_STRING << "dataset '" << fstr.str()
                  << "' can't open. func_name='" << err.getFuncName()
                  << "', msg='" << err.getDetailMsg() << "'." << ENDL;
        return;
    }
    dataset->write(matrix, H5::PredType::STD_I32LE, dataspace3);
    delete dataset;

    // clip-end counts at the detected positions
    fstr.str("/");
    fstr << ref_name << "/ClipEndCount";

    try {
        H5::Exception::dontPrint();
        dataset = new H5::DataSet(file->createDataSet(fstr.str().c_str(),
                H5::PredType::STD_I32LE, dataspace3, ds_creatplist));
    } catch (H5::DataSetIException err) {
        std::cerr << ERROR_STRING << "dataset '" << fstr.str()
                  << "' can't open. func_name='" << err.getFuncName()
                  << "', msg='" << err.getDetailMsg() << "'." << ENDL;
        return;
    }
    dataset->write(clipend_matrix, H5::PredType::STD_I32LE, dataspace3);
    delete dataset;

    delete group;
}
//------------------------------------------------------------------------------
inline void print_usage(const char *cmd) {
    std::cerr << USAGE_STRING << cmd
              << " (-t num_threads=8) -i [bam_fn] -o [matrix_fn]" << ENDL;
}
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
    // parse arguments
    char option;
    std::string input_fn = "", output_fn = "";
    int num_threads = NUM_THREADS;
    while ((option = getopt(argc, argv, "i:o:t:")) != -1) {
        switch (option) {
            case 'i':
                input_fn = optarg;
                break;
            case 'o':
                output_fn = optarg;
                break;
            case 't':
                num_threads = std::atoi(optarg);
                break;
            default:
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
                break;
        }
    }
    if (input_fn.empty() || output_fn.empty()) {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // get RefVector to know about the ref sequences
    BamTools::RefVector refvector = get_refvector(input_fn);
    if (0 >= refvector.size()) {
        std::cerr << ERROR_STRING << "returned invalid refvector." << input_fn
                  << ENDL;
        exit(EXIT_FAILURE);
    }

    // supress output of error messages
    H5::Exception::dontPrint();

    // output file
    H5::H5File *file;
    try {
        H5::Exception::dontPrint();
        file = new H5::H5File(output_fn.c_str(), H5F_ACC_TRUNC);
    } catch (H5::FileIException err) {
        std::cerr << "The output matrix (" << output_fn
                  << ") can't open for writing." << ENDL;
        exit(EXIT_FAILURE);
    }

    // thread args
    pthread_t *thid = new pthread_t[num_threads];
    ThreadCountParam *param = new ThreadCountParam[num_threads];
    IntMatrixType **read_matrix = new IntMatrixType *[num_threads];
    IntMatrixType **clipend_matrix = new IntMatrixType *[num_threads];

    int th_count = 0;
    for (BamTools::RefVector::const_iterator ref = refvector.begin();
            ref != refvector.end(); ++ref) {
        // Counting matrices
        read_matrix[th_count] = new IntMatrixType[ref->RefLength];
        clipend_matrix[th_count] = new IntMatrixType[ref->RefLength];

        // Thread params
        param[th_count].input_fn = input_fn;
        param[th_count].ref_name = ref->RefName;
        param[th_count].ref_start = 0;
        param[th_count].ref_end = ref->RefLength;
        param[th_count].ref_length = ref->RefLength;
        param[th_count].unique_read_count = 0;
        param[th_count].matrix = read_matrix[th_count];
        param[th_count].clipend_matrix = clipend_matrix[th_count];

#pragma omp parallel for
        for (int i = 0; i < ref->RefLength; ++i) {
            read_matrix[th_count][i] = 0;
            clipend_matrix[th_count][i] = 0;
        }

        // counting
        pthread_create(
                &thid[th_count], NULL, thread_make_matrix, &param[th_count]);
        ++th_count;

        if (th_count == num_threads) {
            for (int i = 0; i < th_count; ++i)
                pthread_join(thid[i], NULL);
            for (int i = 0; i < th_count; ++i) {
                write_hdf(file, param[i].ref_name.c_str(), read_matrix[i],
                        clipend_matrix[i], &param[i].ref_length,
                        &param[i].unique_read_count);
                delete[] read_matrix[i];
                delete[] clipend_matrix[i];
            }
            th_count = 0;
        }
    }

    // purge leftovers
    if (0 < th_count) {
        for (int i = 0; i < th_count; ++i)
            pthread_join(thid[i], NULL);
        for (int i = 0; i < th_count; ++i) {
            write_hdf(file, param[i].ref_name.c_str(), read_matrix[i],
                    clipend_matrix[i], &param[i].ref_length,
                    &param[i].unique_read_count);
            delete[] read_matrix[i];
            delete[] clipend_matrix[i];
        }
    }

    delete[] thid;
    delete[] param;
    delete file;
    exit(EXIT_SUCCESS);
}
//------------------------------------------------------------------------------
