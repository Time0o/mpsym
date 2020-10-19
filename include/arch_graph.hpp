#ifndef GUARD_ARCH_GRAPH_H
#define GUARD_ARCH_GRAPH_H

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "arch_graph_system.hpp"
#include "nauty_graph.hpp"
#include "perm_group.hpp"

namespace mpsym
{

namespace internal { class NautyGraph; }

class ArchGraph : public ArchGraphSystem
{
  friend std::ostream &operator<<(std::ostream &os, ArchGraph &ag);

  typedef std::vector<std::string>::size_type processor_type_size_type;
  typedef std::vector<std::string>::size_type channel_type_size_type;

  typedef boost::vecS vertex_selector;
  typedef boost::vecS edge_selector;

  struct VertexProperty { processor_type_size_type type; };
  struct EdgeProperty { channel_type_size_type type; };

  typedef boost::adjacency_list<
    edge_selector, vertex_selector, boost::directedS,
    VertexProperty, EdgeProperty> adjacency_type;

  typedef adjacency_type::vertices_size_type vertices_size_type;
  typedef adjacency_type::edges_size_type edges_size_type;

public:
  typedef processor_type_size_type ProcessorType;
  typedef channel_type_size_type ChannelType;

  typedef std::string const &ProcessorLabel;
  typedef std::string const &ChannelLabel;

  ArchGraph(bool directed = false)
  : _directed(directed)
  {}

  virtual ~ArchGraph() = default;

  std::string to_gap() const override;
  std::string to_json() const override;

  ProcessorType new_processor_type(ProcessorLabel pl = "");
  ChannelType new_channel_type(ChannelLabel cl = "");

  unsigned add_processor(ProcessorType pt);
  unsigned add_processor(ProcessorLabel pl = "");

  unsigned add_processors(unsigned n, ProcessorType pt);
  unsigned add_processors(unsigned n, ProcessorLabel pl = "");

  void add_channel(unsigned pe1, unsigned pe2, ChannelType ct);
  void add_channel(unsigned pe1, unsigned pe2, ChannelLabel cl = "");

  void fully_connect(ChannelType ct);
  void fully_connect(ChannelLabel cl = "");

  bool directed() const { return _directed; }
  unsigned num_processors() const override;
  unsigned num_channels() const override;

private:
  internal::PermGroup automorphisms_(
    AutomorphismOptions const *options) override
  { return automorphisms_nauty(options); }

  bool edge_exists(unsigned from, unsigned to, ChannelType ct) const;
  bool edge_exists_directed(unsigned from, unsigned to, ChannelType ct) const;
  bool edge_exists_undirected(unsigned from, unsigned to, ChannelType ct) const;

  internal::NautyGraph graph_nauty() const;

  std::string to_gap_nauty() const;

  internal::PermGroup automorphisms_nauty(
    AutomorphismOptions const *options) const;

  void dump_processors(std::ostream& os) const;
  void dump_channels(std::ostream& os) const;
  void dump_automorphisms(std::ostream& os);

  adjacency_type _adj;
  bool _directed;

  std::vector<std::string> _processor_types;
  std::vector<std::string> _channel_types;

  std::vector<vertices_size_type> _processor_type_instances;
  std::vector<edges_size_type> _channel_type_instances;
};

std::ostream &operator<<(std::ostream &os, ArchGraph &ag);

}

#endif // GUARD_ARCH_GRAPH_H
