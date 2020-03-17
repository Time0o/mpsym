local mpsym = require 'mpsym'

local processors = {}
local channels = {}

local socs = {}

for i = 1,4 do
  socs[i] = mpsym.append_processors(processors, mpsym.identical_processors(16, 'P'))
end

for i = 1,4 do
  mpsym.append_channels(channels, mpsym.grid_channels(socs[i], 'optical'))
end

for i = 1,3 do
  mpsym.append_channels(channels, mpsym.fully_cross_connected_channels(socs[i], socs[i + 1], 'wireless'))
end


return mpsym.ArchGraph:create{
  processors = processors,
  channels = channels
}
