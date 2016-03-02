#include "parse_region.h"

int parse_region(region_string_map& region_files, const std::string& input_string) {
  std::vector<std::string> pair;
  boost::split(pair, input_string, boost::is_any_of(":"));

  int id = boost::lexical_cast<int>(pair[1]);
  region_files[id] = pair[0];

  return id;
}
