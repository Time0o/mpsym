local mpsym = require 'mpsym'

local args = mpsym.parse_args(args, 'i')

local processors = mpsym.identical_processors(args[1], 'P')
local channels = mpsym.fully_connected_channels(processors, 'C')

return mpsym.ArchGraph:create{
  processors = processors,
  channels = channels
}
