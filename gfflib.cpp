#include "gfflib.h"

//------------------------------------------------------------------------------
bool GffRecord::parse(const std::string &line, bool isParseAttributes = true) {
    hi::StringArray items;
    hi::split(items, line, '\t');
    if (EX_GFFLIB_COL_ATTRIBUTES >= items.size()) {
        std::cerr << ERROR_STRING << "invalid record structure. line=" << line
                  << ENDL;
        return false;
    }

    this->seqid = items.at(EX_GFFLIB_COL_SEQID);
    this->source = items.at(EX_GFFLIB_COL_SOURCE);
    this->type = items.at(EX_GFFLIB_COL_TYPE);
    this->start = std::atoi(items.at(EX_GFFLIB_COL_START).c_str());
    this->end = std::atoi(items.at(EX_GFFLIB_COL_END).c_str());
    this->score = std::atof(items.at(EX_GFFLIB_COL_SCORE).c_str());
    this->strand = items.at(EX_GFFLIB_COL_STRAND)[0];
    this->phase = items.at(EX_GFFLIB_COL_PHASE);
    this->attributes = items.at(EX_GFFLIB_COL_ATTRIBUTES);

    if (isParseAttributes)
        hi::extract_and_add_key_value_pairs(
                attributeItems, this->attributes, ';', '=');

    return true;
}
//------------------------------------------------------------------------------
GffRecord::GffRecord() {
    seqid = "";
}
//------------------------------------------------------------------------------
GffRecord::GffRecord(const std::string &line, bool isParseAttributes = true) {
    this->parse(line, isParseAttributes);
}
//------------------------------------------------------------------------------
bool get_attribute_item(
        const std::string &line, const std::string &key, std::string &value) {
    // do not search by 'key=' to avoid accidental partial match
    // check if key matches the first attribute item
    std::stringstream sstr;
    sstr << key << '=';
    const std::string headPattern = sstr.str();
    if (headPattern == line.substr(0, headPattern.length())) {
        const size_t end = line.find_first_of(";", headPattern.length());
        if (std::string::npos != end) {
            value = line.substr(
                    headPattern.length(), end - headPattern.length());
        } else {
            value = line.substr(headPattern.length());
        }
        return true;
    }

    // second and thereafter attribute item always begins with ';key='
    sstr.str("");
    sstr << ';' << key << '=';
    const std::string keyword = sstr.str();
    const size_t hit = line.find(keyword);
    if (std::string::npos != hit) {
        const size_t end = line.find_first_of(";", hit + keyword.length());
        if (std::string::npos != end) {
            value = line.substr(
                    hit + keyword.length(), end - hit - keyword.length());
        } else {
            value = line.substr(hit + keyword.length());
        }
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
bool replace_attribute_item(
        std::string &line, const std::string &key, const std::string &value) {
    // do not search by 'key=' to avoid accidental partial match
    // check if key matches the first attribute item
    std::stringstream sstr;
    sstr << key << '=';
    const std::string headPattern = sstr.str();
    if (headPattern == line.substr(0, headPattern.length())) {
        const size_t end = line.find_first_of(";", headPattern.length());
        if (std::string::npos != end) {
            line.replace(
                    headPattern.length(), end - headPattern.length(), value);
        } else {
            line.replace(headPattern.length(),
                    line.length() - headPattern.length(), value);
        }
        return true;
    }

    // second and thereafter attribute item always begins with ';key='
    sstr.str("");
    sstr << ';' << key << '=';
    const std::string keyword = sstr.str();
    const size_t hit = line.find(keyword);
    if (std::string::npos != hit) {
        const size_t end = line.find_first_of(";", hit + keyword.length());
        if (std::string::npos != end) {
            line.replace(hit + keyword.length(), end - hit - keyword.length(),
                    value);
        } else {
            line.replace(hit + keyword.length(),
                    line.length() - hit - keyword.length(), value);
        }
        return true;
    }

    // 'key' does not exist -> add one to the end
    std::stringstream attr;
    if (!line.empty() && ';' != line[line.length() - 1]) {
        attr << line << ';';
    }
    attr << key << '=' << value;
    line = attr.str();
    return true;
}
//------------------------------------------------------------------------------
bool remove_attribute_item(std::string &line, const std::string &key) {
    const std::string empty = "";
    // do not search by 'key=' to avoid accidental partial match
    // check if key matches the first attribute item
    std::stringstream sstr;
    sstr << key << '=';
    const std::string headPattern = sstr.str();
    if (headPattern == line.substr(0, headPattern.length())) {
        const size_t end = line.find_first_of(";", headPattern.length());
        if (std::string::npos != end) {
            line.replace(0, end + 1, empty);
        } else {
            line = "";
        }
        return true;
    }

    // second and thereafter attribute item always begins with ';key='
    sstr.str("");
    sstr << ';' << key << '=';
    const std::string keyword = sstr.str();
    const size_t hit = line.find(keyword);
    if (std::string::npos != hit) {
        const size_t end = line.find_first_of(";", hit + keyword.length());
        if (std::string::npos != end) {
            line.replace(hit, end - hit, empty);
        } else {
            line = line.substr(0, hit - 1);
        }
        return true;
    }
    return true;
}
//------------------------------------------------------------------------------
