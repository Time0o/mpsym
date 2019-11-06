#include <iostream>

#include "arch_graph.h"

using cgtl::ArchGraph;

ArchGraph make_haec_plane()
{ return ArchGraph::hyper_mesh(4, 4); }

int main()
{
  ArchGraph haec_plane(make_haec_plane());
}
