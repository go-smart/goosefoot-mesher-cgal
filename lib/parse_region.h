#ifndef PARSE_REGION_H
#define PARSE_REGION_H

typedef std::map< int, std::string > region_string_map;

int parse_region(region_string_map& region_files, const std::string& input_string);

#endif
