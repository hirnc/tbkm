#include "histd.h"

namespace hi {

//------------------------------------------------------------------------------
bool split(StringArray &result, const std::string &line, const char delimiter,
        const SplitMode mode) {
    std::stringstream sstr(line);
    std::string element;
    while (std::getline(sstr, element, delimiter)) {
        if (HISTD_SPLITMODE_ADDALL == mode
                || (HISTD_SPLITMODE_NOEMPTY == mode && !element.empty())) {
            result.push_back(element);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
bool add_record(
        KeyValueDB &targetDB, const std::string &key, const std::string value) {
    KeyValueDB::iterator dbHit = targetDB.find(key);
    if (targetDB.end() != dbHit)
        dbHit->second = value;
    else
        targetDB.insert(std::pair<std::string, std::string>(key, value));
    return true;
}
//------------------------------------------------------------------------------
bool get_record(const KeyValueDB &inputDB, const std::string &key,
        std::string &result) {
    KeyValueDB::const_iterator dbHit = inputDB.find(key);
    if (inputDB.end() != dbHit)
        result = dbHit->second;
    else {
        result = EX_KVDB_ITEM_NULL;
        return false;
    }

    return true;
}
//------------------------------------------------------------------------------
bool add_key_value_pair(
        KeyValueDB &targetDB, const std::string &kvpair, const char delimiter) {
    size_t hitpos = kvpair.find(delimiter);
    if (std::string::npos == hitpos) {
        add_record(targetDB, kvpair, kvpair);
        return true;
    }

    const std::string key = kvpair.substr(0, hitpos);
    const std::string value = kvpair.substr(hitpos + 1);
    add_record(targetDB, key, value);
    return true;
}
//------------------------------------------------------------------------------
bool extract_key_value_pairs(const std::string &line, const char delimiter,
        hi::StringArray &kvpairs) {
    if (!hi::split(kvpairs, line, delimiter))
        return false;
    return true;
}
//------------------------------------------------------------------------------
bool extract_and_add_key_value_pairs(KeyValueDB &targetDB,
        const std::string &line, const char kvp_separator,
        const char kv_delimiter) {
    hi::StringArray kvpairs;
    if (!extract_key_value_pairs(line, kvp_separator, kvpairs))
        return false;

    for (hi::StringArray::const_iterator iter = kvpairs.begin();
            iter != kvpairs.end(); ++iter) {
        add_key_value_pair(targetDB, *iter, kv_delimiter);
    }
    return true;
}
//------------------------------------------------------------------------------
}  // namespace hi
