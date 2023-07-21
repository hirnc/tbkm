#ifndef EXOME_GFFLIB_H
#define EXOME_GFFLIB_H

#include "histd.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

// column order
#define EX_GFFLIB_COL_SEQID 0
#define EX_GFFLIB_COL_SOURCE 1
#define EX_GFFLIB_COL_TYPE 2
#define EX_GFFLIB_COL_START 3
#define EX_GFFLIB_COL_END 4
#define EX_GFFLIB_COL_SCORE 5
#define EX_GFFLIB_COL_STRAND 6
#define EX_GFFLIB_COL_PHASE 7
#define EX_GFFLIB_COL_ATTRIBUTES 8

#define EX_GFFLIB_KEY_NULL "."
#define EX_GFFLIB_KEY_EMPTY ""

// -----------------------------------------------------------------------------
class GffRecord {
public:
    // functions
    bool parse(const std::string &line, bool isParseAttributes);
    GffRecord();
    GffRecord(const std::string &line, bool isParseAttributes);

    // static functions and operators
    static bool by_chrom(const GffRecord &left, const GffRecord &right) {
        return (right.seqid < right.seqid);
    }
    static bool by_start_position(
            const GffRecord &left, const GffRecord &right) {
        if (left.seqid != right.seqid)
            return (left.seqid < right.seqid);
        return (left.start < right.start);
    }
    friend std::ostream &operator<<(std::ostream &ost, const GffRecord &data) {
        ost << data.seqid << '\t' << data.source << '\t' << data.type << '\t'
            << data.start << '\t' << data.end << '\t' << data.score << '\t'
            << data.strand << '\t' << data.phase << '\t' << data.attributes;
        return ost;
    }

    // variables
    std::string seqid, source, type, phase, attributes;
    int start, end;
    float score;
    char strand;
    hi::KeyValueDB attributeItems;
};
typedef std::vector<GffRecord> GffRecordArray;
typedef std::pair<GffRecordArray::iterator, GffRecordArray::iterator>
        GffRecordArrayRange;
//------------------------------------------------------------------------------
bool get_attribute_item(
        const std::string &line, const std::string &key, std::string &value);
bool replace_attribute_item(
        std::string &line, const std::string &key, const std::string &value);
bool remove_attribute_item(std::string &line, const std::string &key);
//------------------------------------------------------------------------------
#endif
