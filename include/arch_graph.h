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
  typedef processor_type_size_type ProcessorType;
  typedef channel_type_size_type ChannelType;

  ProcessorType new_processor_type(std::string label);
  ChannelType new_channel_type(std::string label);

  std::size_t add_processor(ProcessorType pe);
  void add_channel(std::size_t pe1, std::size_t pe2, ChannelType ch);

  std::size_t num_processors() const;
  std::size_t num_channels() const;

  PermGroup generate_automorphisms();
  TaskMapping mapping(std::vector<std::size_t> const &task_assignment) const;
  bool equivalent(TaskMapping const &tm1, TaskMapping const &tm2) const;

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
  friend class ArchGraph;

  std::vector<std::size_t> get() const { return _mapping; }

private:
  TaskMapping(std::vector<std::size_t> mapping, std::size_t eq_class)
    : _mapping(mapping), _eq_class(eq_class) {}

  std::vector<std::size_t> _mapping;
  std::size_t _eq_class;
};

}

#endif // _GUARD_ARCH_GRAPH_H
