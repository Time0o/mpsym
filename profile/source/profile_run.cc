#include <cstdlib>
#include <fstream>
#include <string>

#include <unistd.h>

#include "profile_run.h"
#include "profile_timer.h"
#include "profile_utility.h"

bool run_gap(std::string const &script, double *t)
{
  char ftmp[] = {'X', 'X', 'X', 'X', 'X', 'X'};
  if (mkstemp(ftmp) == -1) {
    error("failed to create temporary file");
    return false;
  }

  class FileRemover
  {
  public:
    FileRemover(std::string const &f) : _f(f) {}
    ~FileRemover() { std::remove(_f.c_str()); }
  private:
    std::string _f;
  } file_remover(ftmp);

  std::ofstream f(ftmp);
  if (f.fail()) {
    error("failed to create temporary file");
    return false;
  }

  f << script;

  f.flush();

  pid_t maybe_child;
  switch ((maybe_child = timer_start())) {
  case -1:
    error("failed to fork child process");
    return false;
    break;
  case 0:
    {
      if (execlp("gap", "gap", "--nointeract", "-q", ftmp, nullptr) == -1) {
        error("failed to exec gap");
        _Exit(EXIT_FAILURE);
      }

      _Exit(EXIT_SUCCESS);
    }
    break;
  }

  *t = timer_stop(maybe_child);

  return true;
}
