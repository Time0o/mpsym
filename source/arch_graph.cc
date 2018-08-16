#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "boost/graph/adjacency_list.hpp"
#include "boost/range/iterator_range_core.hpp"

extern "C" {
  #include "lua.h"
  #include "lualib.h"
  #include "lauxlib.h"
}

#include "arch_graph.h"
#include "dbg.h"
#include "nauty.h"
#include "perm.h"
#include "perm_group.h"

namespace cgtl
{

ArchGraph::ProcessorType ArchGraph::new_processor_type(std::string const &label)
{
  auto id = _processor_types.size();
  _processor_types.push_back(label);
  _processor_type_instances.push_back(0u);

  return id;
}

ArchGraph::ChannelType ArchGraph::new_channel_type(std::string const &label)
{
  auto id = _channel_types.size();
  _channel_types.push_back(label);
  _channel_type_instances.push_back(0u);

  return id;
}

std::size_t ArchGraph::add_processor(ProcessorType pt)
{
  _automorphisms_valid = false;

  _processor_type_instances[pt]++;

  VertexProperty vp {pt};
  return static_cast<std::size_t>(boost::add_vertex(vp, _adj));
}

void ArchGraph::add_channel(std::size_t from, std::size_t to, ChannelType cht)
{
  _automorphisms_valid = false;

  _channel_type_instances[cht]++;

  EdgeProperty ep {cht};
  boost::add_edge(from, to, ep, _adj);
}

std::size_t ArchGraph::num_processors() const
{
  return static_cast<size_t>(boost::num_vertices(_adj));
}

std::size_t ArchGraph::num_channels() const
{
  return static_cast<size_t>(boost::num_edges(_adj));
}

static int generator_degree;
static std::vector<Perm> generators;

static Perm to_perm(int *perm, int degree)
{
  std::vector<unsigned> tmp(degree);
  for (int i = 0; i < degree; ++i)
    tmp[i] = perm[i] + 1;

  return Perm(tmp);
}

static void save_generator(int, int *perm, int *, int, int, int)
{
  generators.push_back(to_perm(perm, generator_degree));
}

static void nauty_free()
{
  naugraph_freedyn();
  nautil_freedyn();
  nauty_freedyn();
}

void ArchGraph::complete()
{
  Dbg(Dbg::DBG) << "=== Determining architecture graph automorphisms";

  // allocate nauty structures
  DYNALLSTAT(graph, g, g_sz);
  DYNALLSTAT(int, lab, lab_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);

  static DEFAULTOPTIONS_GRAPH(options);
  options.defaultptn = FALSE;
  options.userautomproc = save_generator;

  statsblk stats;

  int cts = static_cast<int>(_channel_types.size());
  int cts_log2 = 0; while (cts >>= 1) ++cts_log2;

  int n_orig = static_cast<int>(boost::num_vertices(_adj));
  int n = static_cast<int>(n_orig * (cts_log2 + 1u));
  int m = SETWORDSNEEDED(n);

#ifndef NDBUG
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
#endif

  DYNALLOC2(graph, g, g_sz, m, n, "malloc graph");
  DYNALLOC1(int, lab, lab_sz, n, "malloc lab");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc ptn");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc orbits");

  // construct nauty graph
  Dbg(Dbg::TRACE) << "=== Constructing isomorphic vertex colored graph "
                  << "(" << n << " vertices)";

  EMPTYGRAPH(g, m, n);

  std::vector<vertices_size_type> processor_type_offsets(
    _processor_types.size());

  std::vector<vertices_size_type> processor_type_counters(
    _processor_types.size(), 0u);

  vertices_size_type offs_accu = 0u;
  for (processor_type_size_type i = 0u; i < _processor_types.size(); ++i) {
    processor_type_offsets[i] = offs_accu;
    offs_accu += _processor_type_instances[i];
  }

  /* node numbering:
   *  ...     ...           ...
   *   |       |             |
   * (n+1)---(n+2)-- ... --(n+n)
   *   |       |             |
   *  (1)-----(2)--- ... ---(n)
   */
  for (int level = 0; level <= cts_log2; ++level) {

    if (level > 0) {
      for (int v = 0; v < n_orig; ++v) {
        Dbg(Dbg::TRACE) << "Adding (vertical) edge: " << v + level * n_orig
                        << " => " << v + (level - 1) * n_orig;

        ADDONEEDGE(g, v + level * n_orig, v + (level - 1) * n_orig, m);
      }
    }

    for (auto e : boost::make_iterator_range(boost::edges(_adj))) {
      int t = static_cast<int>(_adj[e].type) + 1;

      if (t & (1 << level)){
        int source = static_cast<int>(boost::source(e, _adj));
        int target = static_cast<int>(boost::target(e, _adj));

        Dbg(Dbg::TRACE) << "Adding (horizontal) edge: "
                        << source + level * n_orig << " => "
                        << target + level * n_orig;

        ADDONEEDGE(g, source + level * n_orig, target + level * n_orig, m);
      }
      t >>= 1u;
    }
  }

  for (int level = 0; level <= cts_log2; ++level) {
    std::fill(processor_type_counters.begin(),
              processor_type_counters.end(), 0u);

    for (int v = 0; v < n_orig; ++v) {
      processor_type_size_type t = _adj[v].type;
      vertices_size_type offs =
        processor_type_offsets[t] + processor_type_counters[t];

      lab[offs + level * n_orig] = v + level * n_orig;

      ptn[offs + level * n_orig] =
        ++processor_type_counters[t] != _processor_type_instances[t];
    }
  }

  Dbg(Dbg::TRACE) << "Colored vertices";

  // call nauty
  generators.clear();
  generator_degree = n_orig;

  densenauty(g, lab, ptn, orbits, &options, &stats, m, n, nullptr);

  DYNFREE(g, g_sz);
  DYNFREE(lab, lab_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  nauty_free();

  if (generators.empty())
    _automorphisms = PermGroup(generator_degree, {Perm(generator_degree)});
  else
    _automorphisms = PermGroup(generator_degree, generators);

  Dbg(Dbg::DBG) << "=== Result";
  Dbg(Dbg::DBG) << _automorphisms;

  _automorphisms_valid = true;
}

PermGroup ArchGraph::automorphisms() const
{
  assert(_automorphisms_valid);

  return _automorphisms;
}

static std::vector<unsigned> min_elem_bruteforce(PermGroup const &ag,
  std::vector<unsigned> const &tasks)
{
  std::vector<unsigned> min_element(tasks);

  Dbg(Dbg::TRACE) << "=== Performing brute force mapping";

  for (Perm const &perm : ag) {
    bool minimal = true;

    for (auto i = 0u; i < tasks.size(); ++i) {
      unsigned permuted = perm[tasks[i] + 1] - 1;

      if (permuted < min_element[i])
        break;

      if (permuted > min_element[i]) {
        minimal = false;
        break;
      }
    }

    if (minimal) {
      for (auto i = 0u; i < tasks.size(); ++i) {
        min_element[i] = perm[tasks[i] + 1] - 1;
      }
    }
  }

  Dbg(Dbg::TRACE) << "Found minimal orbit element: "
                  << tasks << " => " << min_element;

  return min_element;
}

static std::vector<unsigned> min_elem_approx(PermGroup const &ag,
  std::vector<unsigned> const &tasks)
{
  std::vector<Perm> generators(ag.bsgs().sgs());
  std::vector<unsigned> min_element(tasks);

  bool stationary, new_minimum;

  Dbg(Dbg::TRACE) << "=== Performing approximate mapping";

  do {
    stationary = true;

    for (Perm const &gen : generators) {
      new_minimum = false;

      for (unsigned task : min_element) {
        unsigned permuted = gen[task + 1] - 1;

        if (permuted < task) {
          new_minimum = true;
          break;
        }

        if (permuted > task)
          break;
      }

      if (new_minimum) {
        for (unsigned &task : min_element)
          task = gen[task + 1] - 1;

        stationary = false;
        break;
      }
    }
  } while (!stationary);

  Dbg(Dbg::TRACE) << "Found minimal orbit element: "
                  << tasks << " => " << min_element;

  return min_element;
}

TaskMapping ArchGraph::mapping(std::vector<unsigned> const &tasks,
 MappingVariant mapping_variant) const
{
  assert(_automorphisms_valid);
  assert(boost::num_vertices(_adj) > 0u);

  switch (mapping_variant) {
    case MAP_APPROX:
      return TaskMapping(tasks, min_elem_approx(_automorphisms, tasks));
    default:
      return TaskMapping(tasks, min_elem_bruteforce(_automorphisms, tasks));
  }
}

static void lua_parse_err(lua_State *L, std::string const &infile,
  std::string const &err)
{
  Dbg(Dbg::WARN) << "Failed to parse '" << infile << "': " << err;

  lua_close(L);

  throw std::runtime_error("malformed architecture graph description");
}

ArchGraph ArchGraph::fromlua(std::string const &infile)
{
  ArchGraph ag;

  int lua_err;

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  lua_err = luaL_loadfile(L, infile.c_str());
  if (lua_err) {
    Dbg(Dbg::WARN) << "Failed to open '" << infile << "'";

    lua_close(L);

    throw std::runtime_error("failed to open architecture graph description");
  }

  lua_err = lua_pcall(L, 0, 0, 0);
  if (lua_err)
    lua_parse_err(L, infile, lua_tostring(L, -1));

  std::unordered_map<int, std::size_t> pes;
  std::unordered_map<std::string, ProcessorType> pe_types;
  std::unordered_map<std::string, ChannelType> ch_types;

  // parse 'processors' table
  if (lua_getglobal(L, "processors") != LUA_TTABLE)
    lua_parse_err(L, infile, "no 'processors' table defined");

  lua_pushvalue(L, -1);
  lua_pushnil(L);

  while (lua_next(L, -2)) {
    if (!lua_istable(L, -1)) {
      lua_parse_err(L, infile, "malformed element in 'processors' table");

    } else {
      if (!lua_isnumber(L, -2) || !lua_istable(L, -1))
        lua_parse_err(L, infile, "malformed element in 'processors' table");

      lua_pushnil(L);

      if (!lua_next(L, -2) || !lua_isnumber(L, -2) || !lua_isnumber(L, -1))
        lua_parse_err(L, infile, "malformed element in 'processors' table");

      int pe = lua_tonumber(L, -1);
      if (pes.find(pe) != pes.end()) {
        std::stringstream ss;
        ss <<  "processing element " << pe
           << " defined twice in 'processors' table";
        lua_parse_err(L, infile, ss.str());
      }

      lua_pop(L, 1);

      if (!lua_next(L, -2) || !lua_isnumber(L, -2) || !lua_isstring(L, -1))
        lua_parse_err(L, infile, "malformed element in 'processors' table");

      std::string pe_type(lua_tostring(L, -1));
      if (pe_types.find(pe_type) == pe_types.end())
        pe_types[pe_type] = ag.new_processor_type(pe_type);

      pes[pe] = ag.add_processor(pe_types[pe_type]);

      lua_pop(L, 3);
    }
  }

  // parse 'channels' table
  if (lua_getglobal(L, "channels") != LUA_TTABLE)
    lua_parse_err(L, infile, "no 'channels' table defined");

  lua_pushvalue(L, -1);
  lua_pushnil(L);

  while (lua_next(L, -2)) {
    if (!lua_istable(L, -1)) {
      lua_parse_err(L, infile, "malformed element in 'channels' table");

    } else {
      if (!lua_isnumber(L, -2) || !lua_istable(L, -1))
        lua_parse_err(L, infile, "malformed element in 'channels' table");

      lua_pushnil(L);

      if (!lua_next(L, -2) || !lua_isnumber(L, -2) || !lua_isnumber(L, -1))
        lua_parse_err(L, infile, "malformed element in 'channels' table");

      int pe1 = lua_tonumber(L, -1);
      if (pes.find(pe1) == pes.end()) {
        std::stringstream ss;
        ss << "processing element  " << pe1
           << " used in 'channels' table not defined in 'processors' table";
        lua_parse_err(L, infile, ss.str());
      }

      lua_pop(L, 1);

      if (!lua_next(L, -2) || !lua_isnumber(L, -2) || !lua_isnumber(L, -1))
        lua_parse_err(L, infile, "malformed element in 'channels' table");

      int pe2 = lua_tonumber(L, -1);
      if (pes.find(pe2) == pes.end()) {
        std::stringstream ss;
        ss << "processing element " << pe2
           << " used in 'channels' table not defined in 'processors' table";
        lua_parse_err(L, infile, ss.str());
      }

      lua_pop(L, 1);

      if (!lua_next(L, -2) || !lua_isnumber(L, -2) || !lua_isstring(L, -1))
        lua_parse_err(L, infile, "malformed element in 'channels' table");

      std::string ch_type(lua_tostring(L, -1));
      if (ch_types.find(ch_type) == ch_types.end())
        ch_types[ch_type] = ag.new_channel_type(ch_type);

      ag.add_channel(pes[pe1], pes[pe2], ch_types[ch_type]);

      lua_pop(L, 3);
    }
  }

  lua_close(L);

  return ag;
}

void ArchGraph::todot(std::string const &outfile) const
{
  std::ofstream out(outfile);

  if (!out.is_open())
    throw std::runtime_error("failed to create architecture dotfile");

  static char const * const COLORSCHEME = "accent";
  static const unsigned COLORS = 8;
  static char const * const NODESTYLE = "filled";
  static const unsigned LINEWIDTH = 2;

  assert(_processor_types.size() < COLORS
         && "distinguishably many processor types in dot output");

  assert(_channel_types.size() < COLORS
         && "distinguishably many channel types in dot output");

  // construct dotfile...
  out << "graph {\n";
  out << "layout=neato\n";
  out << "splines=true\n";
  out << "overlap=scalexy\n";
  out << "sep=1\n";

  // add vertices
  for (auto v : boost::make_iterator_range(boost::vertices(_adj))) {
    out << v << " [label=" << "PE" << v + 1u << ",style=" << NODESTYLE
        << ",colorscheme=" << COLORSCHEME << COLORS << ",fillcolor="
        << _adj[v].type + 1u << "]\n";
  }

  // add edges
  for (auto e : boost::make_iterator_range(boost::edges(_adj))) {
    auto source = boost::source(e, _adj);
    auto target = boost::target(e, _adj);

    out << source << " -- " << target << " [penwidth=" << LINEWIDTH
        << ",colorscheme=" << COLORSCHEME << COLORS << ",color="
        << _adj[e].type + 1u << "]\n";
  }

  out << "}\n";
}

}
