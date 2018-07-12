#ifndef _GUARD_ARCH_GRAPH_H
#define _GUARD_ARCH_GRAPH_H

#include <cstdlib>
#include <ostream>
#include <string>

#include "boost/graph/adjacency_list.hpp"

#include "perm_group.h"

namespace cgtl
{

struct TaskMapping;

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

  struct VertexProperty { processor_type_size_type type; };
  struct EdgeProperty { channel_type_size_type type; };

  typedef boost::adjacency_list<
    edge_selector, vertex_selector, boost::undirectedS,
    VertexProperty, EdgeProperty> adjacency_type;

public:
  enum MappingVariant { MAP_AUTO, MAP_BRUTEFORCE, MAP_APPROX };

  typedef processor_type_size_type ProcessorType;
  typedef channel_type_size_type ChannelType;

  ProcessorType new_processor_type(std::string label);
  ChannelType new_channel_type(std::string label);

  std::size_t add_processor(ProcessorType pe);
  void add_channel(std::size_t pe1, std::size_t pe2, ChannelType ch);

  std::size_t num_processors() const;
  std::size_t num_channels() const;

  PermGroup generate_automorphisms();
  TaskMapping mapping(std::vector<unsigned> const &tasks,
    MappingVariant mapping_variant = MAP_AUTO) const;

  bool fromlua(std::string const &infile);
  bool todot(std::string const &outfile) const;

private:
  adjacency_type _adj;
  PermGroup _automorphisms;

  std::vector<std::string> _processor_types;
  std::vector<std::string> _channel_types;

  std::vector<vertices_size_type> _processor_type_instances;
  std::vector<edges_size_type> _channel_type_instances;
};

struct TaskMapping {
  TaskMapping(std::vector<unsigned> map, std::vector<unsigned> eq)
    : _mapping(map), _equivalence_class(eq) {}

  std::vector<unsigned> mapping() const { return _mapping; }
  std::vector<unsigned> equivalence_class() const { return _equivalence_class; }

private:
  std::vector<unsigned> _mapping, _equivalence_class;
};

}

#endif // _GUARD_ARCH_GRAPH_H
