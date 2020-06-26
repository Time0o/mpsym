#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

extern "C" {
  #include "lua.h"
  #include "lualib.h"
  #include "lauxlib.h"
}

#include "arch_graph.hpp"
#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"

namespace
{

// stack access

template<typename T,
         typename std::enable_if<std::is_same<std::string, T>::value, int>::type = 0>
T lua_get(lua_State *L, int index = -1)
{
  assert(lua_isstring(L, index));

  return lua_tostring(L, index);
}

template<typename T,
         typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
T lua_get(lua_State *L, int index = -1)
{
  assert(lua_isnumber(L, index));

  lua_Number ret(lua_tonumber(L, index));

  assert(std::floor(ret) == ret
         && ret >= static_cast<lua_Number>(std::numeric_limits<T>::min())
         && ret <= static_cast<lua_Number>(std::numeric_limits<T>::max()));

  return ret;
}

template<typename T,
         typename std::enable_if<std::is_pointer<T>::value, int>::type = 0>
T lua_get(lua_State *L, int index = -1)
{
  assert(lua_isuserdata(L, index));

  return static_cast<T>(lua_touserdata(L, index));
}

template<typename T>
T lua_get_and_pop(lua_State *L)
{
  T ret = lua_get<T>(L);

  lua_pop(L, 1);

  return ret;
}

template<typename T>
T lua_get_from_array(lua_State *L, lua_Integer i, int index = -1)
{
  lua_pushinteger(L, i);

  lua_gettable(L, index - 1);

  return lua_get_and_pop<T>(L);
}

template<typename FUNC>
void lua_foreach_in_table(lua_State *L, FUNC &&f, int index = -1)
{
  lua_pushnil(L);

  while (lua_next(L, index - 1))
    f(L);
}

template<typename FUNC>
void lua_foreach_in_field(lua_State *L,
                          std::string const &field,
                          FUNC &&f,
                          int index = -1)
{
  lua_getfield(L, index, field.c_str());

  lua_foreach_in_table(L, f, -1);

  lua_pop(L, 1);
}

std::string lua_metaname(lua_State *L, int index)
{
  lua_getmetatable(L, index);
  lua_getfield(L, -1, "metaname");

  assert(lua_isstring(L, -1));
  auto ret(lua_get<std::string>(L));

  lua_pop(L, 2);

  return ret;
}

// debugging

// error handling

struct lua_Error : public std::runtime_error
{
  explicit lua_Error(std::string const &what)
  : std::runtime_error("lua: " + what)
  {}
};

// ArchGraphSystem specific functions

bool lua_is_arch_graph(lua_State *L, int index)
{ return lua_metaname(L, index) == "ArchGraph"; }

bool lua_is_arch_graph_cluster(lua_State *L, int index)
{ return lua_metaname(L, index) == "ArchGraphCluster"; }

bool lua_is_arch_uniform_super_graph(lua_State *L, int index)
{ return lua_metaname(L, index) == "ArchUniformSuperGraph"; }

bool lua_is_arch_graph_system(lua_State *L, int index)
{
  return lua_is_arch_graph(L, index)
         || lua_is_arch_graph_cluster(L, index)
         || lua_is_arch_uniform_super_graph(L, index);
}

std::shared_ptr<mpsym::ArchGraphSystem> lua_make_arch_graph_system(lua_State *L);

std::shared_ptr<mpsym::ArchGraph> lua_make_arch_graph(lua_State *L)
{
  assert(lua_is_arch_graph(L, -1));

  auto ag(std::make_shared<mpsym::ArchGraph>());

  // add processors
  std::map<lua_Integer, unsigned> processors;
  std::map<std::string, mpsym::ArchGraph::ProcessorType> processor_types;

  lua_foreach_in_field(L, "_processor_types", [&](lua_State *L){
    auto pl(lua_get_and_pop<std::string>(L));
    processor_types[pl] = ag->new_processor_type(pl);
  });

  lua_foreach_in_field(L, "processors", [&](lua_State *L){
    auto pid(lua_get_from_array<lua_Integer>(L, 1));
    auto pl(lua_get_from_array<std::string>(L, 2));
    lua_pop(L, 1);

    processors[pid] = ag->add_processor(processor_types[pl]);
  });

  // add channels
  std::map<std::string, mpsym::ArchGraph::ChannelType> channel_types;

  lua_foreach_in_field(L, "_channel_types", [&](lua_State *L){
    auto cl(lua_get_and_pop<std::string>(L));
    channel_types[cl] = ag->new_channel_type(cl);
  });

  lua_foreach_in_field(L, "channels", [&](lua_State *L){
    auto source(lua_get_from_array<lua_Integer>(L, 1));
    auto target(lua_get_from_array<lua_Integer>(L, 2));
    auto cl(lua_get_from_array<std::string>(L, 3));
    lua_pop(L, 1);

    ag->add_channel(processors[source], processors[target], channel_types[cl]);
  });

#ifndef NDEBUG
  lua_getfield(L, -1, "_num_processors");
  assert(ag->num_processors() == lua_get_and_pop<unsigned>(L));

  lua_getfield(L, -1, "_num_channels");
  assert(ag->num_channels() == lua_get_and_pop<unsigned>(L));
#endif

  return ag;
}

std::shared_ptr<mpsym::ArchGraphCluster> lua_make_arch_graph_cluster(
  lua_State *L)
{
  auto agc(std::make_shared<mpsym::ArchGraphCluster>());

  lua_foreach_in_table(L, [&](lua_State *L){
    agc->add_subsystem(lua_make_arch_graph_system(L));
  });

  return agc;
}

std::shared_ptr<mpsym::ArchUniformSuperGraph> lua_make_arch_uniform_super_graph(
  lua_State *L)
{
  lua_getfield(L, -1, "super_graph");
  auto super_graph(lua_make_arch_graph_system(L));
  lua_pop(L, 1);

  lua_getfield(L, -1, "proto");
  auto proto(lua_make_arch_graph_system(L));
  lua_pop(L, 1);

  return std::make_shared<mpsym::ArchUniformSuperGraph>(super_graph, proto);
}

std::shared_ptr<mpsym::ArchGraphSystem> lua_make_arch_graph_system(lua_State *L)
{
  assert(lua_is_arch_graph_system(L, -1));

  if (lua_is_arch_graph(L, -1))
    return lua_make_arch_graph(L);
  else if (lua_is_arch_graph_cluster(L, -1))
    return lua_make_arch_graph_cluster(L);
  else if (lua_is_arch_uniform_super_graph(L, -1))
    return lua_make_arch_uniform_super_graph(L);
  else
    throw std::logic_error("unreachable");
}

} // anonymous namespace

namespace mpsym
{

using namespace internal;

std::shared_ptr<ArchGraphSystem> ArchGraphSystem::from_lua_file(
  std::string const &lua_file,
  std::vector<std::string> const &args)
{
  std::ifstream lua_stream(lua_file);

  if (!lua_stream)
    throw std::runtime_error("failed to read lua file");

  std::string lua((std::istreambuf_iterator<char>(lua_stream)),
                   std::istreambuf_iterator<char>());

  return from_lua(lua, args);
}

std::shared_ptr<ArchGraphSystem> ArchGraphSystem::from_lua(
  std::string const &lua,
  std::vector<std::string> const &args)
{
  using namespace std::placeholders;

  // initialize
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  // populate args table
  if (args.size() > 0u) {
    lua_createtable(L, args.size(), 0);

    for (auto i = 0u; i < args.size(); ++i) {
      lua_pushinteger(L, i + 1u);
      lua_pushstring(L, args[i].c_str());
      lua_settable(L, -3);
    }

    lua_setglobal(L, "args");
  }

  // load chunk
  switch (luaL_loadstring(L, lua.c_str())) {
    case LUA_OK:
      break;
    case LUA_ERRSYNTAX:
      throw lua_Error("syntax error while loading chunk");
    case LUA_ERRMEM:
      throw lua_Error("memory error while loading chunk");
    case LUA_ERRGCMM:
      throw lua_Error("garbage collector error while loading chunk");
  }

  // run chunk
  lua_call(L, 0, LUA_MULTRET);

  if (lua_gettop(L) != 1)
    throw lua_Error("chunk did not return singular value");

  // construct ArchGraphSystem
  if (!lua_is_arch_graph_system(L, -1))
    throw lua_Error("invalid ArchGraphSystem descriptor");

  return lua_make_arch_graph_system(L);
}

} // namespace mpsym
