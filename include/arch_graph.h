#ifndef _GUARD_ARCH_GRAPH_H
#define _GUARD_ARCH_GRAPH_H

#include <ostream>
#include <string>

#include "boost/functional/hash.hpp"
#include "boost/graph/adjacency_list.hpp"

#include "perm_group.h"

namespace cgtl
{

class ArchGraph
{
  typedef boost::vecS vertex_selector;
  typedef boost::vecS edge_selector;

  typedef std::vector<std::string>::size_type processor_type_size_type;
  typedef std::vector<std::string>::size_type channel_type_size_type;

  typedef boost::adjacency_list<
    edge_selector, vertex_selector, boost::bidirectionalS> traits;

  typedef traits::vertices_size_type vertices_size_type;
  typedef traits::edges_size_type edges_size_type;

  struct VertexProperty {
    processor_type_size_type type;
  };

  struct EdgeProperty {
    channel_type_size_type type;
  };

  typedef boost::adjacency_list<
    edge_selector, vertex_selector, boost::undirectedS,
    VertexProperty, EdgeProperty> adjacency_type;

public:
  enum Autom { AUTOM_TOTAL, AUTOM_PROCESSORS, AUTOM_CHANNELS };

  typedef processor_type_size_type ProcessorType;
  typedef channel_type_size_type ChannelType;

  typedef adjacency_type::vertex_descriptor Processor;
  typedef adjacency_type::edge_descriptor Channel;
  typedef std::vector<Processor> Mapping;

  struct MappingRepresentation {
    MappingRepresentation(Mapping const &m) : _hash(boost::hash_value(m)) {}

    bool operator==(MappingRepresentation const &rhs) {
      return (rhs._hash == _hash);
    }

    bool operator!=(MappingRepresentation const &rhs) {
      return !((*this) == rhs);
    }

  private:
    size_t _hash;
  };

  ProcessorType new_processor_type(std::string label);
  ChannelType new_channel_type(std::string label);

  Processor add_processor(ProcessorType pe);
  void add_channel(ProcessorType pe1, ProcessorType pe2, ChannelType ch);

  size_t num_processors() const {
    return static_cast<size_t>(boost::num_vertices(_adj));
  }
  size_t num_channels() const {
    return static_cast<size_t>(boost::num_edges(_adj));
  }

  PermGroup automorphisms(Autom at = AUTOM_TOTAL) const;

  MappingRepresentation mapping_representation(Mapping const &m) const;
  bool mappings_equivalent(Mapping const &m1, Mapping const &m2) const;
  bool mappings_equivalent(MappingRepresentation const &m1,
                           MappingRepresentation const &m2) const;

  bool fromlua(std::string const &infile);
  bool todot(std::string const &outfile) const;

private:
  PermGroup processor_automorphisms() const;
  PermGroup channel_automorphisms() const;

  adjacency_type _adj;

  std::vector<std::string> _processor_types;
  std::vector<std::string> _channel_types;

  std::vector<vertices_size_type> _processor_type_instances;
  std::vector<edges_size_type> _channel_type_instances;
};

}

#endif // _GUARD_ARCH_GRAPH
