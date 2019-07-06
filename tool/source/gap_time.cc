#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <utility>

#include <libgen.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/wait.h>
#include <unistd.h>


enum { CHAR_GAP_READY = 'g', CHAR_GAP_START = 's' };

namespace {
  template<typename Arg, typename... Args>
  void error(Arg &&arg, Args&&... args) {
    std::cerr << std::forward<Arg>(arg);
    (void)(int[]){(std::cerr << ' ' << std::forward<Args>(args), 0)...};
  }
}

int main(int argc, char **argv) {
  if (argc != 2) {
    error("usage:", basename(argv[0]), "GAP_SCRIPT");
    exit(EXIT_FAILURE);
  }

  pid_t child;
  switch((child = fork())) {
  case -1:
    error("failed to fork child process");
    exit(EXIT_FAILURE);
    break;
  case 0:
    {
      if (execlp("gap", "gap", "--nointeract", "-q", argv[1], nullptr) == -1) {
        error("failed to exec gap");
        _Exit(EXIT_FAILURE);
      }
    }
    break;
  default:
    {
      int status;
      waitpid(child, &status, 0);

      if (!WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS)
        exit(EXIT_FAILURE);

      struct tms tms;
      times(&tms);

      double t = static_cast<double>(tms.tms_cutime) / sysconf(_SC_CLK_TCK);

      std::cout << "cpu time: " << t << "s\n";
    }
  }

  exit(EXIT_SUCCESS);
}
