local mpsym = require 'mpsym'

local processors = {}
local channels = {}

local socs = {}

for i = 1,16 do
  socs[i] = mpsym.append_processors(processors, mpsym.identical_processors(16, 'P'))
end

for i = 1,16 do
  mpsym.append_channels(channels, mpsym.fully_connected_channels(socs[i], 'shared memory'))
end

for i = 1,4 do
  for j = 1,4 do
    local k = (i - 1) * 4 + j

    if j < 4 then
      mpsym.append_channels(channels, mpsym.fully_cross_connected_channels(socs[k], socs[k + 1], 'C'))
    end

    if i < 4 then
      mpsym.append_channels(channels, mpsym.fully_cross_connected_channels(socs[k], socs[k + 4], 'C'))
    end
  end
end

return mpsym.ArchGraph:create{
  processors = processors,
  channels = channels
}
