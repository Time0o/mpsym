#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <string>

#include "gmock/gmock.h"

#include "dbg.hpp"

#define USAGE "USAGE: TEST [OPTIONS]\n" \
              " -h           display this message\n" \
              " -o TESTCASE  only run the testcase TESTCASE\n" \
              " -v           increase verbosity, can be passed multiple times\n"

int main(int argc, char** argv) {
  int verbosity = 0;

  option const longopts[] {
    {"help", no_argument, nullptr, 'h'},
    {"only", required_argument, nullptr, 'o'}
  };

  int opt, index;
  while ((opt = getopt_long(argc, argv, "ho:v", longopts, &index)) != -1) {
    switch(opt) {
      case 'h':
        std::cout << USAGE;
        exit(EXIT_SUCCESS);
      case 'o':
        testing::GTEST_FLAG(filter) = std::string("*") + optarg + "*";
        break;
      case 'v':
        ++verbosity;
        break;
      case '?':
        std::cerr << USAGE;
        exit(EXIT_FAILURE);
    }
  }

  switch(verbosity) {
    case 0:
      DBG_SET_LOGLEVEL(WARN);
      break;
    case 1:
      DBG_SET_LOGLEVEL(INFO);
      break;
    case 2:
      DBG_SET_LOGLEVEL(DEBUG);
      break;
    default:
      DBG_SET_LOGLEVEL(TRACE);
      break;
  }

  testing::InitGoogleMock(&argc, argv);
  return RUN_ALL_TESTS();
}
