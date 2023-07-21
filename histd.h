#ifndef HISTD_H
#define HISTD_H

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <new>
#include <sstream>
#include <string>
#include <vector>

#ifdef NONUNIX
#include <ctype>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

// cerr strings
#define ERROR_STRING "\x1b[31;1;7mError: \x1b[m"
#define WARNING_STRING "\x1b[32;1mWarning: \x1b[m"
#define INFO_STRING "\x1b[36;1mInfo: \x1b[m"
#define USAGE_STRING "\x1b[36;1mUsage: \x1b[m"

#define ENDL '\n'
#define EX_KVDB_ITEM_NULL "."
//-----------------------------------------------------------------------------
namespace hi {
typedef std::vector<std::string> StringArray;
typedef std::map<std::string, std::string> KeyValueDB;


enum SplitMode { HISTD_SPLITMODE_ADDALL, HISTD_SPLITMODE_NOEMPTY };
bool split(StringArray &result, const std::string &line, const char delimiter,
        const SplitMode mode = HISTD_SPLITMODE_ADDALL);
bool add_record(
        KeyValueDB &targetDB, const std::string &key, const std::string value);
bool get_record(
        const KeyValueDB &inputDB, const std::string &key, std::string &result);
bool add_key_value_pair(
        KeyValueDB &targetDB, const std::string &kvpair, const char delimiter);
bool extract_key_value_pairs(const std::string &line, const char delimiter,
        hi::StringArray &kvpairs);
bool extract_and_add_key_value_pairs(KeyValueDB &targetDB,
        const std::string &line, const char kvp_separator,
        const char kv_delimiter);
}  // namespace hi

#endif
