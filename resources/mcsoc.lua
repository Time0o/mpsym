processors = {
  {0, 'A7'},
  {1, 'A7'},
  {2, 'A7'},
  {3, 'A7'},
  {4, 'A15'},
  {5, 'A15'},
  {6, 'A15'},
  {7, 'A15'}
}

channels = {}

for i=0,7 do
  table.insert(channels, {i, i, 'L1'})
end

for i=0,3 do
  for j=i,3 do
    table.insert(channels, {i, j, 'L2'})
  end
end

for i=4,7 do
  for j=i,7 do
    table.insert(channels, {i, j, 'L2'})
  end
end

for i=0,7 do
  for j=i,7 do
    table.insert(channels, {i, j, 'SRAM'})
  end
end
