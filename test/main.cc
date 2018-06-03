#include <cstring>
#include <iostream>

#include "dbg.h"
#include "gmock/gmock.h"

static void usage_err(char const *progname)
{
  std::cerr << "USAGE: " << progname << " [-v|-vv|-vvv]\n";
}

int main(int argc, char** argv) {
  bool run = true;

  if (argc > 2) {
     usage_err(argv[0]);
     run = false;
  } else if (argc == 2) {
    if (strcmp(argv[1], "-v") == 0) {
      Dbg::loglevel = Dbg::INFO;
    } else if (strcmp(argv[1], "-vv") == 0) {
      Dbg::loglevel = Dbg::DBG;
    } else if (strcmp(argv[1], "-vvv") == 0) {
      Dbg::loglevel = Dbg::TRACE;
    } else {
      usage_err(argv[0]);
      run = false;
    }
  }

  if (run) {
    ::testing::InitGoogleMock(&argc, argv);
    return RUN_ALL_TESTS();
  }

  return -1;
}
