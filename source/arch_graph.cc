#include <algorithm>
#include <cassert>
#include <fstream>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>

extern "C" {
  #include "lua.h"
  #include "lualib.h"
  #include "lauxlib.h"
}

#include "arch_graph.h"
#include "dbg.h"
#include "nauty.h"
#include "partial_perm.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm.h"
#include "perm_group.h"

namespace cgtl
{

namespace
{

int generator_degree;
std::vector<Perm> generators;

Perm to_perm(int *perm, int degree)
{
  std::vector<unsigned> tmp(degree);
  for (int i = 0; i < degree; ++i)
    tmp[i] = perm[i] + 1;

  return Perm(tmp);
}

void save_generator(int, int *perm, int *, int, int, int)
{
  generators.push_back(to_perm(perm, generator_degree));
}

void nauty_free()
{
  naugraph_freedyn();
  nautil_freedyn();
  nauty_freedyn();
}

std::vector<unsigned> min_elem_bruteforce(
  PermGroup const &ag, std::vector<unsigned> const &tasks,
  unsigned min_pe, unsigned max_pe)
{
  std::vector<unsigned> min_element(tasks);

  Dbg(Dbg::DBG) << "Performing brute force mapping";

  for (Perm const &perm : ag) {
    bool minimal = true;

    for (auto i = 0u; i < tasks.size(); ++i) {
      unsigned task_base = tasks[i];
      if (task_base < min_pe || task_base > max_pe)
        continue;

      unsigned task = task_base - min_pe;

      unsigned permuted = perm[task + 1u] - 1u + min_pe;

      if (permuted < min_element[i])
        break;

      if (permuted > min_element[i]) {
        minimal = false;
        break;
      }
    }

    if (minimal) {
      for (auto i = 0u; i < tasks.size(); ++i) {
        unsigned task_base = tasks[i];
        if (task_base < min_pe || task_base > max_pe)
          continue;

        unsigned task = task_base - min_pe;

        min_element[i] = perm[task + 1u] - 1u + min_pe;
      }
    }
  }

  Dbg(Dbg::DBG) << "Found minimal orbit element: " << min_element;

  return min_element;
}

std::vector<unsigned> min_elem_approx(
  PermGroup const &ag, std::vector<unsigned> const &tasks,
  unsigned min_pe, unsigned max_pe)
{
  std::vector<Perm> generators(ag.bsgs().strong_generators);
  std::vector<unsigned> min_element(tasks);

  bool stationary, new_minimum;

  Dbg(Dbg::TRACE) << "Performing approximate mapping";

  do {
    stationary = true;

    for (Perm const &gen : generators) {
      new_minimum = false;

      for (unsigned task_base : min_element) {
        if (task_base < min_pe || task_base > max_pe)
          continue;

        unsigned task = task_base - min_pe;

        unsigned permuted = gen[task + 1u] - 1u;

        if (permuted < task) {
          new_minimum = true;
          break;
        }

        if (permuted > task)
          break;
      }

      if (new_minimum) {
        for (unsigned &task_base : min_element) {
          if (task_base < min_pe || task_base > max_pe)
            continue;

          unsigned task = task_base - min_pe;

          task_base = gen[task + 1u] - 1u + min_pe;
        }

        stationary = false;
        break;
      }
    }
  } while (!stationary);

  Dbg(Dbg::DBG) << "Found minimal orbit element: " << min_element;

  return min_element;
}

} // anonymous namespace

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

unsigned ArchGraph::add_processor(ProcessorType pt)
{
  _automorphisms_valid = false;

  _processor_type_instances[pt]++;

  VertexProperty vp {pt};
  return static_cast<unsigned>(boost::add_vertex(vp, _adj));
}

void ArchGraph::add_channel(unsigned from, unsigned to, ChannelType cht)
{
  _automorphisms_valid = false;

  _channel_type_instances[cht]++;

  EdgeProperty ep {cht};
  boost::add_edge(from, to, ep, _adj);
}

unsigned ArchGraph::num_processors() const
{
  return static_cast<unsigned>(boost::num_vertices(_adj));
}

unsigned ArchGraph::num_channels() const
{
  return static_cast<size_t>(boost::num_edges(_adj));
}

void ArchGraph::complete()
{
  if (_automorphisms_valid)
    return;

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

PartialPermInverseSemigroup ArchGraph::partial_automorphisms() const
{
  unsigned n = num_processors();

  std::function<bool(PartialPerm const &)>
  is_partial_automorphism = [&](PartialPerm const &pperm) {
    for (unsigned x : pperm.dom()) {
      unsigned y = pperm[x];
      if (_adj[x - 1u].type != _adj[y - 1u].type)
        return false;

      for (unsigned z : pperm.dom()) {
        auto edge(boost::edge(x - 1u, z - 1u, _adj));
        if (!std::get<1>(edge))
          continue;

        auto edge_perm(boost::edge(y - 1u, z - 1u, _adj));
        if (!std::get<1>(edge_perm))
          return false;

        auto edge_type = _adj[std::get<0>(edge)].type;
        auto edge_perm_type = _adj[std::get<0>(edge_perm)].type;

        if (edge_type != edge_perm_type)
          return false;
      }
    }

    return true;
  };

  struct Domain {
    Domain(std::vector<bool> const &set, unsigned limit)
      : set(set), limit(limit) {}

    std::vector<bool> set;
    unsigned limit;
  };

  PartialPermInverseSemigroup res;

  std::function<void(Domain const &, std::vector<unsigned> const &)>
  backtrack = [&](Domain const &domain, std::vector<unsigned> const &image) {
    std::vector<unsigned> tmp(domain.limit, 0u);

    unsigned j = 0u;
    for (unsigned i = 0u; i < domain.limit; ++i) {
      if (domain.set[i])
        tmp[i] = image[j++];
    }

    PartialPerm pperm(tmp);
    Dbg(Dbg::TRACE) << "Considering " << pperm << ':';

    if (is_partial_automorphism(pperm)) {
      Dbg(Dbg::TRACE) << "=> Is a partial automorphism";
      if (!pperm.id()) {
        Dbg(Dbg::TRACE) << "==> Adjoining new generator";
        res.adjoin({pperm}, true);
      }
    } else {
      Dbg(Dbg::TRACE)
        << "=> Is not a partial automorphism, pruning backtrack tree...";
      return;
    }

    Domain next_domain(domain);

    std::vector<unsigned> next_image(image);
    next_image.resize(next_image.size() + 1u);

    for (unsigned i = domain.limit; i < n; ++i) {
      next_domain.set[i] = true;
      next_domain.limit = i + 1u;

      for (unsigned j = 1u; j <= n; ++j) {
        if (std::find(image.begin(), image.end(), j) != image.end())
          continue;

        next_image.back() = j;
        backtrack(next_domain, next_image);
      }

      next_domain.set[i] = false;
    }
  };

  Dbg(Dbg::DBG)
    << "Finding partial automorphisms for arch graph with automorphism group:";
  Dbg(Dbg::DBG) << _automorphisms;

  backtrack(Domain(std::vector<bool>(num_processors(), false), 0u), {});

  return res;
}

TaskMapping ArchGraph::mapping(
 std::vector<unsigned> const &tasks, unsigned offset,
 MappingVariant mapping_variant) const
{
  Dbg(Dbg::DBG) << "Requested task mapping for: " << tasks;

  unsigned min_pe = offset;
  unsigned max_pe = offset + num_processors() - 1u;

  assert(_automorphisms_valid);
  assert(boost::num_vertices(_adj) > 0u);

#ifndef NDEBUG
  if (min_pe != 0u) {
    Dbg(Dbg::TRACE) << "Mapping shifted range ["
                    << min_pe << ", " << max_pe << "]";
  }
#endif

  switch (mapping_variant) {
    case MAP_APPROX:
      return TaskMapping(tasks, min_elem_approx(_automorphisms, tasks,
                                                min_pe, max_pe));
    default:
      return TaskMapping(tasks, min_elem_bruteforce(_automorphisms, tasks,
                                                    min_pe, max_pe));
  }
}

/*
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

  std::unordered_map<int, unsigned> pes;
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
*/

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

void ArchGraphCluster::add_subsystem(
  std::shared_ptr<ArchGraphSystem> const &ags)
{
  _subsystems.push_back(ags);
}

unsigned ArchGraphCluster::num_processors() const
{
  unsigned res = 0u;
  for (auto const &subsystem : _subsystems)
    res += subsystem->num_processors();

  return res;
}

unsigned ArchGraphCluster::num_channels() const
{
  unsigned res = 0u;
  for (auto const &subsystem : _subsystems)
    res += subsystem->num_channels();

  return res;
}

void ArchGraphCluster::complete()
{
  assert(!_subsystems.empty());

  for (auto const &subsystem : _subsystems)
    subsystem->complete();

  _automorphisms = _subsystems[0]->automorphisms();
  for (auto i = 1u; i < _subsystems.size(); ++i) {
    _automorphisms = PermGroup::direct_product(
      _automorphisms, _subsystems[i]->automorphisms());
  }

  _automorphisms_valid = true;
}

PermGroup ArchGraphCluster::automorphisms() const
{
  assert(_automorphisms_valid);

  return _automorphisms;
}

TaskMapping ArchGraphCluster::mapping(
  std::vector<unsigned> const &tasks, unsigned offset,
  MappingVariant mapping_variant) const
{
  Dbg(Dbg::DBG) << "Requested task mapping for: " << tasks;

  assert(_subsystems.size() > 0u);

  TaskMapping res(tasks, tasks);

  unsigned offs = offset;
  for (auto i = 0u; i < _subsystems.size(); ++i) {
    unsigned next_offs = offs + _subsystems[i]->num_processors();

    Dbg(Dbg::DBG) << "Subsystem (no. " << i << ", "
                  << "pe's " << offs << "-" << next_offs - 1u << ")";

    res = _subsystems[i]->mapping(res.equivalence_class(), offs, mapping_variant);

    Dbg(Dbg::DBG) << "Yields: " << res.equivalence_class();

    offs = next_offs;
  }

  return TaskMapping(tasks, res.equivalence_class());
}

} // namespace cgtl
