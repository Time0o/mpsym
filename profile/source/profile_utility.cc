#include <string>

#include "profile_utility.h"

#define RESOURCE_DIR "../resources/"

std::string resource_path(std::string const &resource)
{
  return RESOURCE_DIR + resource;
}
