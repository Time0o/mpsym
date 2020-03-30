local mpsym = require 'mpsym'

local processors = {}
local p1_processors = mpsym.append_processors(processors, mpsym.identical_processors(4, 'P1'))
local p2_processors = mpsym.append_processors(processors, mpsym.identical_processors(4, 'P2'))

local channels = {}
--mpsym.append_channels(channels, mpsym.self_connected_channels(processors, 'L1'))
--mpsym.append_channels(channels, mpsym.self_connected_channels(p1_processors, 'L2'))
--mpsym.append_channels(channels, mpsym.self_connected_channels(p2_processors, 'L2'))
--mpsym.append_channels(channels, mpsym.fully_connected_channels(processors, 'RAM'))
mpsym.append_channels(channels, mpsym.self_connected_channels(processors, 'C'))
mpsym.append_channels(channels, mpsym.self_connected_channels(p1_processors, 'C'))
mpsym.append_channels(channels, mpsym.self_connected_channels(p2_processors, 'C'))
mpsym.append_channels(channels, mpsym.fully_connected_channels(processors, 'C'))

return mpsym.ArchGraph:create{
  processors = processors,
  channels = channels
}
