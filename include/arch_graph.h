#ifndef _GUARD_ARCH_GRAPH_H
#define _GUARD_ARCH_GRAPH_H

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "arch_graph_system.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"

namespace cgtl
{

using ProcessorLabel = char const *;
using ChannelLabel = char const *;
using SubsystemLabel = char const *;

#define DEFAULT_PROCESSOR_LABEL ""
#define DEFAULT_CHANNEL_LABEL ""
#define DEFAULT_SUBSYSTEM_LABEL ""

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
    edge_selector, vertex_selector, boost::undirectedS,
    VertexProperty, EdgeProperty> adjacency_type;

  typedef adjacency_type::vertices_size_type vertices_size_type;
  typedef adjacency_type::edges_size_type edges_size_type;

public:
  typedef processor_type_size_type ProcessorType;
  typedef channel_type_size_type ChannelType;

  static ArchGraph fully_connected(unsigned n,
                                   ProcessorLabel pl = DEFAULT_PROCESSOR_LABEL,
                                   ChannelLabel cl = DEFAULT_CHANNEL_LABEL);

  static ArchGraph regular_mesh(unsigned width,
                                unsigned height,
                                ProcessorLabel pl = DEFAULT_PROCESSOR_LABEL,
                                ChannelLabel cl = DEFAULT_CHANNEL_LABEL);

  static ArchGraph hyper_mesh(unsigned width,
                              unsigned height,
                              ProcessorLabel pl = DEFAULT_PROCESSOR_LABEL,
                              ChannelLabel cl = DEFAULT_CHANNEL_LABEL);

  ProcessorType new_processor_type(ProcessorLabel pl = DEFAULT_PROCESSOR_LABEL);
  ChannelType new_channel_type(ChannelLabel cl = DEFAULT_CHANNEL_LABEL);

  unsigned add_processor(ProcessorType pe);
  void add_channel(unsigned pe1, unsigned pe2, ChannelType ch);

  unsigned num_processors() const override;
  unsigned num_channels() const override;

  PartialPermInverseSemigroup partial_automorphisms() override;

private:
  PermGroup update_automorphisms() override;

  void create_mesh(unsigned width,
                   unsigned height,
                   ProcessorType pe,
                   ChannelType ch);

  void dump_processors(std::ostream& os) const;
  void dump_channels(std::ostream& os) const;
  void dump_automorphisms(std::ostream& os);

  adjacency_type _adj;

  std::vector<std::string> _processor_types;
  std::vector<std::string> _channel_types;

  std::vector<vertices_size_type> _processor_type_instances;
  std::vector<edges_size_type> _channel_type_instances;
};

std::ostream &operator<<(std::ostream &os, ArchGraph &ag);

}

#endif // _GUARD_ARCH_GRAPH_H
