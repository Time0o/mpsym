#include <algorithm>
#include <cassert>
#include <fstream>
#include <functional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "arch_graph.h"
#include "dbg.h"
#include "nauty.h"
#include "partial_perm.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"

namespace cgtl
{

namespace
{

int nauty_generator_degree;
PermSet nauty_generators;

void nauty_save_generator(int, int *perm, int *, int, int, int)
{
  std::vector<unsigned> tmp(nauty_generator_degree);
  for (int i = 0; i < nauty_generator_degree; ++i)
    tmp[i] = perm[i] + 1;

  nauty_generators.emplace(tmp);
}

void nauty_free()
{
  naugraph_freedyn();
  nautil_freedyn();
  nauty_freedyn();
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
  assert(!_automorphisms_valid);

  _processor_type_instances[pt]++;

  VertexProperty vp {pt};
  return static_cast<unsigned>(boost::add_vertex(vp, _adj));
}

void ArchGraph::add_channel(unsigned from, unsigned to, ChannelType cht)
{
  assert(!_automorphisms_valid);

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
  return static_cast<unsigned>(boost::num_edges(_adj));
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
  options.userautomproc = nauty_save_generator;

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
  nauty_generators.clear();
  nauty_generator_degree = n_orig;

  densenauty(g, lab, ptn, orbits, &options, &stats, m, n, nullptr);

  DYNFREE(g, g_sz);
  DYNFREE(lab, lab_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);
  nauty_free();

  if (nauty_generators.empty())
    _automorphisms = PermGroup(nauty_generator_degree,
                               {Perm(nauty_generator_degree)});
  else
    _automorphisms = PermGroup(nauty_generator_degree, nauty_generators);

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
      return TaskMapping(tasks, min_elem_approx(tasks, min_pe, max_pe));
    default:
      return TaskMapping(tasks, min_elem_bruteforce(tasks, min_pe, max_pe));
  }
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


ArchGraph ArchGraph::fully_connected(
  unsigned n, std::string const &pe_label, std::string const &ch_label)
{
  assert(n > 0u);

  ArchGraph ag;

  auto pe(ag.new_processor_type(pe_label));
  auto ch(ag.new_channel_type(ch_label));

  std::vector<ArchGraph::ProcessorType> processors;
  for (unsigned i = 0u; i < n; ++i)
    processors.push_back(ag.add_processor(pe));

  for (unsigned i = 0u; i < n; ++i) {
    for (unsigned j = 0u; j < n; ++j) {
      if (j == i)
        continue;

      ag.add_channel(processors[i], processors[j], ch);
    }
  }

  ag._automorphisms = PermGroup::symmetric(n);
  ag._automorphisms_valid = true;

  return ag;
}

ArchGraph ArchGraph::regular_mesh(
  unsigned width, unsigned height,
  std::string const &pe_label, std::string const &ch_label)
{
  ArchGraph ag;

  auto pe(ag.new_processor_type(pe_label));
  auto ch(ag.new_channel_type(ch_label));

  ag.create_mesh(width, height, pe, ch);

  if (height == width) {
    ag._automorphisms = PermGroup::dihedral(8);
    ag._automorphisms_valid = true;
  } else {
    ag.complete();
  }

  return ag;
}

ArchGraph ArchGraph::hyper_mesh(
  unsigned width, unsigned height,
  std::string const &pe_label, std::string const &ch_label)
{
  ArchGraph ag;

  auto pe(ag.new_processor_type(pe_label));
  auto ch(ag.new_channel_type(ch_label));

  ag.create_mesh(width, height, pe, ch);

  for (unsigned r = 0u; r < height; ++r) {
    unsigned pe1 = r * width;
    unsigned pe2 = pe1 + width - 1;

    ag.add_channel(pe1, pe2, ch);
  }

  for (unsigned c = 0u; c < width; ++c) {
    unsigned pe1 = c;
    unsigned pe2 = pe1 + (height - 1) * width;

    ag.add_channel(pe1, pe2, ch);
  }

  ag.complete();

  return ag;
}

std::vector<unsigned> ArchGraph::min_elem_bruteforce(
  std::vector<unsigned> const &tasks,
  unsigned min_pe,
  unsigned max_pe) const
{
  std::vector<unsigned> min_element(tasks);

  Dbg(Dbg::DBG) << "Performing brute force mapping";

  for (Perm const &perm : _automorphisms) {
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

std::vector<unsigned> ArchGraph::min_elem_approx(
  std::vector<unsigned> const &tasks,
  unsigned min_pe,
  unsigned max_pe) const
{
  PermSet generators(_automorphisms.bsgs().strong_generators());
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

void ArchGraph::create_mesh(unsigned width,
                            unsigned height,
                            ProcessorType pe,
                            ChannelType ch)
{
  assert(width > 0u && height > 0u);

  std::vector<std::vector<ArchGraph::ProcessorType>> processors(height);

  for (unsigned i = 0u; i < height; ++i) {
    for (unsigned j = 0u; j < width; ++j) {
      processors[i].push_back(add_processor(pe));
    }
  }

  std::pair<int, int> offsets[] = {{-1, 0}, {0, 1}, {1, 0}, {0, -1}};

  int iheight = static_cast<int>(height);
  int iwidth = static_cast<int>(width);

  for (int i = 0u; i < iheight; ++i) {
    for (int j = 0u; j < iwidth; ++j) {
      for (auto const &offs : offsets) {
        int i_offs = i + offs.first;
        int j_offs = j + offs.second;
        if (i_offs < 0 || i_offs >= iheight || j_offs < 0 || j_offs >= iwidth)
          continue;

        add_channel(processors[i][j], processors[i_offs][j_offs], ch);
      }
    }
  }
}

void ArchGraph::dump_processors(std::ostream& os) const
{
  std::vector<std::vector<unsigned>> pes_by_type(_processor_types.size());

  for (auto pe : boost::make_iterator_range(boost::vertices(_adj))) {
    auto pt = _adj[pe].type;

    pes_by_type[pt].push_back(pe);
  }

  os << "processors: [";

  for (auto pt = 0u; pt < pes_by_type.size(); ++pt) {
    os << "\n  type " << pt;

    auto pt_str = _processor_types[pt];
    if (!pt_str.empty())
      os << " (" << pt_str << ")";

    os << ": [";

    for (auto i = 0u; i < pes_by_type[pt].size(); ++i) {
      os << pes_by_type[pt][i];
      if (i < pes_by_type[pt].size() - 1u)
        os << ", ";
    }

    os << "]";
  }

  os << "\n]";
}

void ArchGraph::dump_channels(std::ostream& os) const
{
  std::vector<std::vector<std::set<unsigned>>> chs_by_type(_channel_types.size());

  for (auto &chs : chs_by_type)
    chs.resize(num_processors());

  for (auto pe1 : boost::make_iterator_range(boost::vertices(_adj))) {
    for (auto e : boost::make_iterator_range(boost::out_edges(pe1, _adj))) {
      auto pe2 = boost::target(e, _adj);
      auto ch = _adj[e].type;

      chs_by_type[ch][pe1].insert(pe2);
    }
  }

  os << "channels: [";

  for (auto ct = 0u; ct < chs_by_type.size(); ++ct) {
    os << "\n  type " << ct;

    auto ct_str = _channel_types[ct];
    if (!ct_str.empty())
      os << " (" << ct_str << ")";

    os << ": [";

    for (auto pe = 0u; pe < chs_by_type[ct].size(); ++pe) {
      auto adj(chs_by_type[ct][pe]);

      if (adj.empty())
        continue;

      os << "\n    " << pe << ": [";

      auto it = adj.begin();
      while (it != adj.end()) {
        os << *it;
        if (++it != adj.end())
          os << ", ";
      }

      os << "]";
    }

    os << "\n  ]";
  }

  os << "\n]";
}

void ArchGraph::dump_automorphisms(std::ostream& os) const
{
  os << "automorphism group: [";

  if (_automorphisms_valid) {
    auto gens(_automorphisms.bsgs().strong_generators());

    for (auto i = 0u; i < gens.size(); ++i) {
      os << "\n  " << gens[i];

      if (i < gens.size() - 1u)
        os << ",";
    }
  } else {
    os << "...";
  }

  os << "\n]";
}

std::ostream& operator<<(std::ostream& os, ArchGraph const &ag)
{
  if (ag.num_processors() == 0u) {
    os << "empty architecture graph";
    return os;
  }

  ag.dump_processors(os);
  os << "\n";
  ag.dump_channels(os);
  os << "\n";
  ag.dump_automorphisms(os);
  os << "\n";

  return os;
}

} // namespace cgtl
