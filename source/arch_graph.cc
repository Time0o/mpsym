#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <cassert>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#include "arch_graph.h"
#include "dbg.h"
#include "nauty.h"
#include "perm.h"
#include "perm_group.h"

namespace cgtl
{

static Perm to_perm(int *perm, int degree)
{
  std::vector<unsigned> tmp(degree);
  for (int i = 0; i < degree; ++i)
    tmp[i] = perm[i] + 1;

  return Perm(tmp);
}

static int generator_degree;
static std::vector<Perm> generators;

static void save_generator(int, int *perm, int *, int, int, int)
{
  generators.push_back(to_perm(perm, generator_degree));
}

// free all dynamic memory allocated by a call to nauty
static void nauty_free()
{
  naugraph_freedyn();
  nautil_freedyn();
  nauty_freedyn();
}

ArchGraph::ProcessorType ArchGraph::new_processor_type(std::string label)
{
  auto id = _processor_types.size();
  _processor_types.push_back(label);
  _processor_type_instances.push_back(0u);

  return id;
}

ArchGraph::ChannelType ArchGraph::new_channel_type(std::string label)
{
  auto id = _channel_types.size();
  _channel_types.push_back(label);
  _channel_type_instances.push_back(0u);

  return id;
}

ArchGraph::Processor ArchGraph::add_processor(ArchGraph::ProcessorType pt)
{
  _processor_type_instances[pt]++;

  VertexProperty vp {pt};
  return boost::add_vertex(vp, _adj);
}

void ArchGraph::add_channel(ArchGraph::Processor from, ArchGraph::Processor to,
                            ArchGraph::ChannelType cht)
{
  _channel_type_instances[cht]++;

  EdgeProperty ep {cht};
  boost::add_edge(from, to, ep, _adj);
}


PermGroup ArchGraph::processor_automorphisms() const
{
  Dbg(Dbg::DBG) << "=== Determining processor automorphisms";

  // allocate nauty structures
  DYNALLSTAT(graph, g, g_sz);
  DYNALLSTAT(int, lab, lab_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);

  static DEFAULTOPTIONS_GRAPH(options);
  options.defaultptn = FALSE;
  options.userautomproc = save_generator;

  statsblk stats;

  size_t n = static_cast<size_t>(boost::num_vertices(_adj));
  size_t m = SETWORDSNEEDED(n);

#ifndef NDBUG
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
#endif

  DYNALLOC2(graph, g, g_sz, m, n, "malloc graph");
  DYNALLOC1(int, lab, lab_sz, n, "malloc lab");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc ptn");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc orbits");

  // construct nauty graph
  Dbg(Dbg::TRACE) << "=== Constructing nauty graph (" << n << " vertices)";

  EMPTYGRAPH(g, m, n);

  for (auto e : boost::make_iterator_range(boost::edges(_adj))) {
    int source = static_cast<int>(boost::source(e, _adj));
    int target = static_cast<int>(boost::target(e, _adj));

    Dbg(Dbg::TRACE) << "Adding edge: " << source << " => " << target;
    ADDONEEDGE(g, source, target, m);
  }

  // color vertices
  std::vector<vertices_size_type> processor_type_offsets(
    _processor_types.size(), 0u);
  std::vector<vertices_size_type> processor_type_counters(
    _processor_types.size(), 0u);

  vertices_size_type offs_accu = 0u;
  for (processor_type_size_type i = 0u; i < _processor_types.size(); ++i) {
    processor_type_offsets[i] = offs_accu;
    offs_accu += _processor_type_instances[i];
  }

  for (auto v : boost::make_iterator_range(boost::vertices(_adj))) {
    processor_type_size_type t = _adj[v].type;
    vertices_size_type offs =
      processor_type_offsets[t] + processor_type_counters[t];

    lab[offs] = static_cast<int>(v);
    ptn[offs] = ++processor_type_counters[t] != _processor_type_instances[t];
  }

  // call nauty
  generators.clear();
  generator_degree = n;

  densenauty(g, lab, ptn, orbits, &options, &stats, m, n, nullptr);

  DYNFREE(g, g_sz);
  DYNFREE(lab, lab_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  nauty_free();

  if (generators.empty()) {
    PermGroup ret(n, {Perm(n)});

    Dbg(Dbg::DBG) << "=== Result";
    Dbg(Dbg::DBG) << "Trivial automorphism group";
    return ret;

  } else {
    PermGroup ret(n, generators);

    Dbg(Dbg::DBG) << "=== Result";
    Dbg(Dbg::DBG) << ret;
    return ret;
  }
}

PermGroup ArchGraph::channel_automorphisms() const
{
  Dbg(Dbg::DBG) << "=== Determining channel automorphisms";

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

  int n_orig = static_cast<int>(boost::num_edges(_adj));
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
  Dbg(Dbg::TRACE) << "=== Constructing isomorphic vertex colored graph (" << n
                  << " vertices)";

  EMPTYGRAPH(g, m, n);

  /* node numbering:
   *  ...     ...           ...
   *   |       |             |
   * (n+1)---(n+2)-- ... --(n+n)
   *   |       |             |
   *  (1)-----(2)--- ... ---(n)
   */
  for (int level = 0; level <= cts_log2; ++level) {

    if (level > 0) {
      for (int i = 0; i < n_orig; ++i) {
        Dbg(Dbg::TRACE) << "Adding (vertical) edge: " << i + level * n_orig
                        << " => " << i + (level - 1) * n_orig;

        ADDONEEDGE(g, i + level * n_orig, i + (level - 1) * n_orig, m);
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

  for (int i = 0; i < n; ++i) {
    lab[i] = i;
    ptn[i] = (i + 1) % n_orig != 0;
  }
  Dbg(Dbg::TRACE) << "Colored vertices";

  generators.clear();
  generator_degree = n_orig;
  densenauty(g, lab, ptn, orbits, &options, &stats, m, n, nullptr);

  DYNFREE(lab, lab_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  DYNFREE(g, g_sz);
  nauty_free();

  if (generators.empty()) {
    PermGroup ret(n_orig, {Perm(n_orig)});

    Dbg(Dbg::DBG) << "=== Result";
    Dbg(Dbg::DBG) << "Trivial automorphism group";
    return ret;

  } else {
    PermGroup ret(n_orig, generators);

    Dbg(Dbg::DBG) << "=== Result";
    Dbg(Dbg::DBG) << ret;
    return ret;
  }
}

PermGroup ArchGraph::automorphisms(ArchGraph::Autom at) const
{
  switch (at) {
    case AUTOM_PROCESSORS:
      return processor_automorphisms();
      break;
    case AUTOM_CHANNELS:
      return channel_automorphisms();
      break;
    default:
      throw std::runtime_error("TODO");
  }
}

void ArchGraph::todot(std::string const &outfile) const
{
  std::ofstream out(outfile);

  if (!out.is_open()) {
    Dbg(Dbg::WARN) << "Could not produce dotfile, failed to open file '"
                   << outfile << "'";
    return;
  }

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
