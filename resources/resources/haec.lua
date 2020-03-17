local mpsym = require 'mpsym'

local super_graph_clusters = mpsym.identical_clusters(4, 'SoC')
local super_graph_channels = mpsym.linear_channels(super_graph_clusters, 'wireless')

local proto_processors = mpsym.identical_processors(16, 'P')
local proto_channels = mpsym.grid_channels(proto_processors, 'optical')

return mpsym.ArchUniformSuperGraph:create{
  super_graph = mpsym.ArchGraph:create{
    clusters = super_graph_clusters,
    channels = super_graph_channels
  },
  proto = mpsym.ArchGraph:create{
    processors = proto_processors,
    channels = proto_channels
  }
}
