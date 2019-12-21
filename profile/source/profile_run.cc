#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include "profile_run.h"
#include "profile_timer.h"

namespace
{
  void dup_fd(int from, int to)
  {
    for (;;) {
      if (dup2(from, to) == -1) {
        if (errno == EINTR)
          continue;
        else
          throw std::runtime_error("dup failed");
      } else {
        break;
      }
    }
  }

  std::string read_output(int from, bool echo)
  {
    static char buf[256];

    std::string res;

    for (;;) {
      auto count = read(from, buf, sizeof(buf));

      if (count == -1) {
        if (errno == EINTR)
          continue;
        else
          throw std::runtime_error("read failed");
        break;
      } else if (count == 0) {
        break;
      } else {
        std::string block(buf, buf + count);

        if (echo)
          std::cout << block << std::flush;

        res += block;
      }
    }

    return res;
  }
}

std::string run_gap(std::string const &script,
                    bool hide_output,
                    bool hide_errors,
                    double *t)
{
  // create temporary gap script

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

  // create pipe for capturing gaps output

  int fds[2];
  if (pipe(fds) == -1)
    throw std::runtime_error("failed to create pipe");

  // run gap in a child process

  timer_start();

  pid_t child;
  switch ((child = fork())) {
  case -1:
    throw std::runtime_error("failed to fork child process");
  case 0:
    {
      dup_fd(fds[1], STDOUT_FILENO);

      fcntl(STDOUT_FILENO, F_SETFL, fcntl(STDOUT_FILENO, F_GETFL) | O_DIRECT);

      close(fds[1]);
      close(fds[0]);

      if (hide_errors) {
        int dev_null = open("/dev/null", O_WRONLY);
        if (dev_null != -1)
          dup_fd(dev_null, STDERR_FILENO);
      }

      if (execlp("gap", "gap", "--nointeract", "-q", ftmp, nullptr) == -1)
        _Exit(EXIT_FAILURE);

      _Exit(EXIT_SUCCESS);
    }
  }

  close(fds[1]);

  auto output(read_output(fds[0], !hide_output));

  close(fds[0]);

  if (t)
    *t = timer_stop(child);

  return output;
}
