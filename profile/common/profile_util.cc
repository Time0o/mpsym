#include <cstring>
#include <ostream>
#include <sstream>
#include <vector>

#include "profile_util.h"
#include "timer.h"

void debug_timer_dump(char const *timer)
{
  std::ostream *os = TIMER_GET_OUT();

  std::stringstream ss;
  TIMER_SET_OUT(&ss);

  TIMER_DUMP(timer);
  debug(ss.str());

  TIMER_SET_OUT(os);
}

std::vector<std::string> split(std::string const &str, char const *delim)
{
  std::vector<std::string> res;

  std::size_t pos = 0u, nextpos;

  for (;;) {
    nextpos = str.find(delim, pos);

    if (nextpos == std::string::npos) {
      res.push_back(str.substr(pos, str.size() - pos));
      break;
    } else {
      res.push_back(str.substr(pos, nextpos - pos));
    }

    pos = nextpos + strlen(delim);
  }

  return res;
}

std::string join(std::vector<std::string> const &strs, char const *delim)
{
  if (strs.empty())
    return "";

  std::string res(strs[0]);
  for (auto i = 1u; i < strs.size(); ++i)
    res += delim + strs[i];

  return res;
}
