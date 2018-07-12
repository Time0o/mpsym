#include <string>

#include "profile_utility.h"

#define RESOURCE_DIR "resources/profile/"

std::string resource_path(std::string const &resource)
{
  return RESOURCE_DIR + resource;
}
