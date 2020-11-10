#ifndef GUARD_ARCH_GRAPH_H
#define GUARD_ARCH_GRAPH_H

#include <cassert>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "arch_graph_system.hpp"
#include "iterator.hpp"
#include "nauty_graph.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

namespace mpsym
{

namespace internal { class NautyGraph; }

class ArchGraph : public ArchGraphSystem
{
  using processor_type_size_type = std::vector<std::string>::size_type;
  using channel_type_size_type = std::vector<std::string>::size_type;

  template<typename L>
  using typed_channel_dict_type = std::unordered_map<
    unsigned,
    std::vector<std::pair<unsigned, L>>
  >;

  using untyped_channel_dict_type = std::unordered_map<
    unsigned,
    std::vector<unsigned>
  >;

  using vertex_selector = boost::vecS;
  using edge_selector = boost::vecS;

  struct VertexProperty { processor_type_size_type type; };
  struct EdgeProperty { channel_type_size_type type; };

  using adjacency_type = boost::adjacency_list<
    edge_selector, vertex_selector, boost::directedS,
    VertexProperty, EdgeProperty>;

  using vertices_size_type = adjacency_type::vertices_size_type;
  using edges_size_type = adjacency_type::edges_size_type;

public:
  using ProcessorType = processor_type_size_type;
  using ChannelType = channel_type_size_type;

  template<typename L = void>
  using ChannelDict = typename std::conditional<
    std::is_void<L>::value,
    untyped_channel_dict_type,
    typed_channel_dict_type<L>
  >::type;

  ArchGraph(bool directed = false)
  : _directed(directed)
  {}

  virtual ~ArchGraph() = default;

  std::string to_gap() const override;
  std::string to_json() const override;

  ProcessorType new_processor_type(std::string const &pl);
  ChannelType new_channel_type(std::string const &cl);

  unsigned add_processor(ProcessorType pt);
  unsigned add_processor(std::string const &pl);

  template<typename T>
  unsigned add_processors(unsigned n, T pt)
  {
    assert(n > 0u);

    for (unsigned i = 0u; i < n - 1u; ++i)
      add_processor(pt);

    return add_processor(pt);
  }

  void add_channel(unsigned pe1, unsigned pe2, ChannelType ct);
  void add_channel(unsigned pe1, unsigned pe2, std::string const &cl);

  template<typename T>
  void add_channels(ChannelDict<T> const &channels)
  {
    for (auto const &tmp : channels) {
      unsigned pe1 = tmp.first;
      auto const &to = tmp.second;

      for (auto const &c : to) {
        unsigned pe2 = c.first;
        T ct = c.second;

        add_channel(pe1, pe2, ct);
      }
    }
  }

  template<typename T>
  void add_channels(ChannelDict<> const &channels, T ct)
  {
    for (auto const &tmp : channels) {
      unsigned pe1 = tmp.first;
      auto const &to = tmp.second;

      for (auto &pe2 : to)
        add_channel(pe1, pe2, ct);
    }
  }

  template<typename T>
  void fully_connect(T ct)
  {
    for (unsigned pe1 = 0u; pe1 < num_processors(); ++pe1) {
      for (unsigned pe2 = (directed() ? 0u : pe1 + 1u); pe2 < num_processors(); ++pe2) {
        if (pe1 != pe2)
          add_channel(pe1, pe2, ct);
      }
    }
  }

  void fully_connect(ProcessorType pt, ChannelType ct);
  void fully_connect(std::string const &pl, std::string const &cl);

  template<typename T>
  void self_connect(T ct)
  {
    for (unsigned pe = 0u; pe < num_processors(); ++pe)
      add_channel(pe, pe, ct);
  }

  void self_connect(ProcessorType pt, ChannelType ct);
  void self_connect(std::string const &pl, std::string const &cl);

  bool directed() const;
  bool effectively_directed() const;

  unsigned num_processors() const override;
  unsigned num_channels() const override;

private:
  internal::PermGroup automorphisms_(
    AutomorphismOptions const *options,
    internal::timeout::flag aborted) override
  { return automorphisms_nauty(options, aborted); }

  void init_repr_(AutomorphismOptions const *options,
                  internal::timeout::flag aborted) override
  { automorphisms(options, aborted); }

  // Convenience functions

  ChannelType assert_channel_type(std::string const &cl);
  ProcessorType assert_processor_type(std::string const &cl);

  bool channel_exists(unsigned from, unsigned to, ChannelType ct) const;
  bool channel_exists_directed(unsigned from, unsigned to, ChannelType ct) const;
  bool channel_exists_undirected(unsigned from, unsigned to, ChannelType ct) const;

  using pe_it = adjacency_type::vertex_iterator;
  using pe = pe_it::value_type;

  boost::iterator_range<pe_it> processors() const
  { return boost::make_iterator_range(boost::vertices(_adj)); }

  processor_type_size_type num_processor_types() const
  { return _processor_types.size(); }

  ProcessorType processor_type(pe pe) const
  { return _adj[pe].type; }

  std::string processor_type_str(pe pe) const
  { return _processor_types[processor_type(pe)]; }

  using ch_it = adjacency_type::edge_iterator;
  using ch_out_it = adjacency_type::out_edge_iterator;
  using ch = ch_it::value_type;

  boost::iterator_range<ch_it> channels() const
  { return boost::make_iterator_range(boost::edges(_adj)); }

  pe source(ch ch) const
  { return boost::source(ch, _adj); }

  pe target(ch ch) const
  { return boost::target(ch, _adj); }

  boost::iterator_range<ch_out_it> out_channels(pe pe) const
  { return boost::make_iterator_range(boost::out_edges(pe, _adj)); }

  channel_type_size_type num_channel_types() const
  { return _channel_types.size(); }

  ChannelType channel_type(ch ch) const
  { return _adj[ch].type; }

  std::string channel_type_str(ch ch) const
  { return _channel_types[channel_type(ch)]; }

  // Nauty

  internal::NautyGraph graph_nauty() const;

  std::string to_gap_nauty() const;

  internal::PermSet automorphism_generators_nauty();

  internal::PermGroup automorphisms_nauty(
    AutomorphismOptions const *options,
    internal::timeout::flag aborted);

  adjacency_type _adj;
  bool _directed;

  std::vector<std::string> _processor_types;
  std::vector<std::string> _channel_types;

  std::vector<vertices_size_type> _processor_type_instances;
  std::vector<edges_size_type> _channel_type_instances;
};

}

#endif // GUARD_ARCH_GRAPH_H
