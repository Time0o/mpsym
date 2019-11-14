#ifndef _GUARD_ARCH_GRAPH_H
#define _GUARD_ARCH_GRAPH_H

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "partial_perm.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_mapping.h"

namespace cgtl
{

class ArchGraphSystem
{
public:
  enum MappingVariant { MAP_AUTO, MAP_BRUTEFORCE, MAP_APPROX };

  virtual unsigned num_processors() const = 0;
  virtual unsigned num_channels() const = 0;

  virtual void complete() = 0;
  virtual bool completed() const = 0;

  virtual PermGroup automorphisms() const = 0;

  virtual PartialPermInverseSemigroup partial_automorphisms() const = 0;

  virtual TaskMapping mapping(
    std::vector<unsigned> const &tasks, unsigned offset = 0u,
    MappingVariant mapping_variant = MAP_AUTO) const = 0;
};

class ArchGraph : public ArchGraphSystem
{
  friend std::ostream &operator<<(std::ostream &os, ArchGraph const &ag);

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

  ProcessorType new_processor_type(std::string const &label = "");
  ChannelType new_channel_type(std::string const &label = "");

  unsigned add_processor(ProcessorType pe);
  void add_channel(unsigned pe1, unsigned pe2, ChannelType ch);

  unsigned num_processors() const override;
  unsigned num_channels() const override;

  void complete() override;
  bool completed() const override { return _automorphisms_valid; }

  PermGroup automorphisms() const override;

  PartialPermInverseSemigroup partial_automorphisms() const override;

  TaskMapping mapping(
    std::vector<unsigned> const &tasks, unsigned offset = 0u,
    MappingVariant mapping_variant = MAP_AUTO) const override;

  static ArchGraph fully_connected(unsigned n,
    std::string const &pe_label = "", std::string const &ch_label = "");

  static ArchGraph regular_mesh(unsigned width, unsigned height,
    std::string const &pe_label = "", std::string const &ch_label = "");

  static ArchGraph hyper_mesh(unsigned width, unsigned height,
    std::string const &pe_label = "", std::string const &ch_label = "");

private:
  std::vector<unsigned> min_elem_bruteforce(std::vector<unsigned> const &tasks,
                                            unsigned min_pe,
                                            unsigned max_pe) const;

  std::vector<unsigned> min_elem_approx(std::vector<unsigned> const &tasks,
                                        unsigned min_pe,
                                        unsigned max_pe) const;

  void create_mesh(unsigned width,
                   unsigned height,
                   ProcessorType pe,
                   ChannelType ch);

  void dump_processors(std::ostream& os) const;
  void dump_channels(std::ostream& os) const;
  void dump_automorphisms(std::ostream& os) const;

  adjacency_type _adj;
  PermGroup _automorphisms;
  bool _automorphisms_valid = false;

  std::vector<std::string> _processor_types;
  std::vector<std::string> _channel_types;

  std::vector<vertices_size_type> _processor_type_instances;
  std::vector<edges_size_type> _channel_type_instances;
};

class ArchGraphCluster : public ArchGraphSystem
{
public:
  void add_subsystem(std::shared_ptr<ArchGraphSystem> const &ags);

  unsigned num_processors() const override;
  unsigned num_channels() const override;

  void complete() override;
  bool completed() const override { return _automorphisms_valid; }

  PermGroup automorphisms() const override;

  PartialPermInverseSemigroup partial_automorphisms() const override {
    throw std::logic_error("not implemented");
  }

  TaskMapping mapping(
    std::vector<unsigned> const &tasks, unsigned offset = 0u,
    MappingVariant mapping_variant = MAP_AUTO) const override;

private:
  PermGroup _automorphisms;
  bool _automorphisms_valid = false;

  std::vector<std::shared_ptr<ArchGraphSystem>> _subsystems;
};

class ArchUniformSuperGraph : public ArchGraphSystem
{
public:
  typedef ArchGraph::ProcessorType SubsystemType;
  typedef ArchGraph::ChannelType SubsystemChannelType;

  ArchUniformSuperGraph() {};

  ArchUniformSuperGraph(
    std::shared_ptr<ArchGraphSystem> const &subsystem_proto,
    std::string const &subsystem_label = "");

  void set_subsystem(
    std::shared_ptr<ArchGraph> const &subsystem_proto,
    std::string const &subsystem_label = "");

  SubsystemChannelType new_subsystem_channel_type(std::string const &label = "");

  SubsystemType add_subsystem();

  void add_subsystem_channel(
    SubsystemType ss1, SubsystemType ss2, SubsystemChannelType ch);

  unsigned num_processors() const override;
  unsigned num_channels() const override;

  void complete() override;
  bool completed() const override { return _automorphisms_valid; }

  PermGroup automorphisms() const override;

  PartialPermInverseSemigroup partial_automorphisms() const override {
    throw std::logic_error("not implemented");
  }

  TaskMapping mapping(
    std::vector<unsigned> const &tasks, unsigned offset = 0u,
    MappingVariant mapping_variant = MAP_AUTO) const override {

    (void)tasks;
    (void)offset;
    (void)mapping_variant;

    throw std::logic_error("not implemented");
  }

private:
  ArchGraph _subsystem_supergraph;
  std::shared_ptr<ArchGraphSystem> _subsystem_proto;
  ArchGraph::ProcessorType _subsystem_processor_type;

  PermGroup _automorphisms;
  bool _automorphisms_valid = false;
};

std::ostream &operator<<(std::ostream &os, ArchGraph const &ag);

}

#endif // _GUARD_ARCH_GRAPH_H
