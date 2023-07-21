#include "histd.h"

#include "gfflib.h"
#include "hdf_base_depth_reader.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#define EX_GFFC_MIN_DEPTH 5
//------------------------------------------------------------------------------
bool set_chromosome(
        HdfBaseDepthReader *hdf, const int nFiles, const char *chromName) {
    for (int i = 0; i < nFiles; ++i) {
        if (!hdf[i].set_target_chromosome(chromName)
                || !hdf[i].set_target_dataset(
                        "BaseDepth", H5::PredType::STD_I32LE))
            return false;
    }
    return true;
}
//------------------------------------------------------------------------------
bool open_hdfs(const hi::StringArray &inputFiles, HdfBaseDepthReader *hdfs) {
    bool isError = false;
    for (hi::StringArray::const_iterator fn = inputFiles.begin();
            fn != inputFiles.end(); ++fn) {
        const int pos = std::distance(inputFiles.begin(), fn);
        if (!hdfs[pos].open(fn->c_str())) {
            std::cerr << WARNING_STRING << "failed to open an input HDF ("
                      << *fn << "). [code: " << __LINE__ << "]" << ENDL;
            isError = true;
        }
    }
    if (isError) {
        return false;
    }
    return true;
}
//------------------------------------------------------------------------------
void init_buffer(int *coveredBases, float *avgDepth, const int nElements,
        const int value = 0) {
    for (int i = 0; i < nElements; ++i) {
        coveredBases[i] = value;
        avgDepth[i] = (float)value;
    }
    return;
}
//------------------------------------------------------------------------------
bool get_cover_stat(HdfBaseDepthReader &hdf, const int *start,
        const int *szRegion, int *coveredBases, const int minDepth,
        float *avgDepth) {
    IntType *matrix = new IntType[*szRegion];
    if (!hdf.get_matrix(start, szRegion, matrix)) {
        std::cerr << WARNING_STRING
                  << "failed to fetch a matrix. start=" << *start
                  << ", size=" << szRegion << ENDL;
        return false;
    }
    // counting
    int totalDP = 0;
    for (int pos = 0; pos < *szRegion; ++pos) {
        if (minDepth <= matrix[pos]) {
            *coveredBases += 1;
        }
        totalDP += matrix[pos];
    }
    *avgDepth = (float)totalDP / (float)*szRegion;
    delete[] matrix;
    return true;
}
//------------------------------------------------------------------------------
bool get_cover_stat(HdfBaseDepthReader *hdfs, const int nFiles,
        const int *start, const int *szRegion, int *coveredBases,
        const int minDepth, float *avgDepth) {
    bool isError = false;
    for (int i = 0; i < nFiles; ++i) {
        if (!get_cover_stat(hdfs[i], start, szRegion, &coveredBases[i],
                    minDepth, &avgDepth[i])) {
            isError = true;
        }
    }
    if (isError) {
        return false;
    }
    return true;
}
//------------------------------------------------------------------------------
bool determine_gff_coverage(HdfBaseDepthReader *hdfs, const int nFiles,
        const GffRecordArray &records, const int minDepth, std::ostream &ofs) {
    const char sep = '\t';
    // common buffer
    int *coveredBases = new int[nFiles];
    float *avgDepth = new float[nFiles];

    // process each record in GFF
    std::string lastChrom = "";
    for (GffRecordArray::const_iterator record = records.begin();
            record != records.end(); ++record) {
        // new chromosome
        if (record->seqid != lastChrom) {
            if (!set_chromosome(hdfs, nFiles, record->seqid.c_str())) {
                std::cerr << WARNING_STRING << "seqid (" << record->seqid
                          << ") does not exist in HDF matrix. Skipped." << ENDL;
                continue;
            }
            lastChrom = record->seqid;
        }
        const int szRegion = record->end - record->start + 1;
        init_buffer(coveredBases, avgDepth, nFiles, 0);
        get_cover_stat(hdfs, nFiles, &(record->start), &szRegion, coveredBases,
                minDepth, avgDepth);
        // write results
        ofs << *record;
        for (int i = 0; i < nFiles; ++i) {
            ofs << sep << avgDepth[i] << sep << coveredBases[i] << sep
                << (float)coveredBases[i] / (float)szRegion;
        }
        ofs << ENDL;
    }
    delete[] coveredBases;
    delete[] avgDepth;
    return true;
}
//------------------------------------------------------------------------------
bool read_gff_from_file(const char *gffFn, GffRecordArray &records) {
    std::ifstream infile(gffFn, std::ios::in);
    if (infile.fail()) {
        std::cerr << ERROR_STRING << "failed to open the input GFF (" << gffFn
                  << "). Aborted. [code: " << __LINE__ << "]" << ENDL;
        return false;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if ('#' == line[0] || line.empty()) {
            continue;
        }
        GffRecord gff(line, false);
        records.push_back(gff);
    }
    infile.close();
    return true;
}
//------------------------------------------------------------------------------
void print_usage(const char *cmd) {
    std::cerr << USAGE_STRING << cmd << " (options) matrix1 matrix2..." << ENDL;
    std::cerr << "Available options:" << ENDL;
    std::cerr << " -i  input GFF filename [MANDATORY]" << ENDL;
    std::cerr << " -m  min read depth to consider 'covered' ["
              << EX_GFFC_MIN_DEPTH << "]" << ENDL << ENDL;
    return;
}
//------------------------------------------------------------------------------
int main(int argc, char **argv) {
    std::string gffFn = "";
    int minDepth = EX_GFFC_MIN_DEPTH;
    // parse arguments
    char option;
    while ((option = getopt(argc, argv, "i:m:h")) != -1) {
        switch (option) {
            case 'i':
                gffFn = optarg;
                break;
            case 'm':
                minDepth = std::atoi(optarg);
                if (0 > minDepth) {
                    std::cerr << WARNING_STRING
                              << "minDepth must bea positive integer or zero. "
                                 "Using a default setting (-m "
                              << EX_GFFC_MIN_DEPTH << ")." << ENDL;
                    minDepth = EX_GFFC_MIN_DEPTH;
                }
                break;
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_SUCCESS);
                break;
            default:
                std::cerr << WARNING_STRING
                          << "unknown option specified and ignored." << ENDL;
                break;
        }
    }

    // check mandatory arguments
    if (gffFn.empty()) {
        std::cerr << ERROR_STRING << "an input GFF (-i) is mandatory." << ENDL
                  << ENDL;
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // create input file list
    hi::StringArray inputFiles, sampleNames;
    for (int i = optind; i < argc; ++i) {
        const char *inputFn = argv[i];
        inputFiles.push_back(inputFn);
        // sample names
        sampleNames.push_back(inputFn);
    }

    // read feature coordinates from GFF
    GffRecordArray records;
    if (!read_gff_from_file(gffFn.c_str(), records)) {
        exit(EXIT_FAILURE);
    }

    // input HDFs
    HdfBaseDepthReader *hdfs = new HdfBaseDepthReader[inputFiles.size()];
    if (!open_hdfs(inputFiles, hdfs)) {
        exit(EXIT_FAILURE);
    }

    // write header
    std::cout << "#CHROM\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattri"
                 "butes";
    for (hi::StringArray::const_iterator name = sampleNames.begin();
            name != sampleNames.end(); ++name) {
        std::cout << '\t' << *name << ".avgDepth\t" << *name
                  << ".coveredBases\t" << *name << ".coveredFrac";
    }
    std::cout << ENDL;

    // process matrices
    if (!determine_gff_coverage(
                hdfs, inputFiles.size(), records, minDepth, std::cout)) {
        exit(EXIT_FAILURE);
    }

    // clean-up
    for (size_t i = 0; i < inputFiles.size(); ++i) {
        hdfs[i].close();
    }
    delete[] hdfs;
    exit(EXIT_SUCCESS);
}
//------------------------------------------------------------------------------
