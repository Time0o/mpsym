#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include "profile_run.h"
#include "profile_timer.h"

void run_gap(std::string const &script, double *t)
{
  char ftmp[] = {'X', 'X', 'X', 'X', 'X', 'X'};
  if (mkstemp(ftmp) == -1)
    throw std::runtime_error("failed to create temporary file");

  class FileRemover
  {
  public:
    FileRemover(std::string const &f) : _f(f) {}
    ~FileRemover() { std::remove(_f.c_str()); }
  private:
    std::string _f;
  } file_remover(ftmp);

  std::ofstream f(ftmp);
  if (f.fail())
    throw std::runtime_error("failed to create temporary file");

  f << script;

  f.flush();

  pid_t maybe_child;
  switch ((maybe_child = timer_start())) {
  case -1:
    throw std::runtime_error("failed to fork child process");
  case 0:
    {
      int dev_null = open("/dev/null", O_WRONLY);
      if (dev_null != -1)
        dup2(dev_null, STDERR_FILENO);

      if (execlp("gap", "gap", "--nointeract", "-q", ftmp, nullptr) == -1)
        _Exit(EXIT_FAILURE);

      _Exit(EXIT_SUCCESS);
    }
    break;
  }

  *t = timer_stop(maybe_child);
}
