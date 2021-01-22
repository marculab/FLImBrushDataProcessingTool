function [map_out] =pad_zero_and_stitch(map_in)
%this function take in a map that has different number of points per line
%and pad the shoet lines with zeros and convert it to an 2D array with
%equal number of lines
% map_in: input map which is a n by 1 cell.
% map_out: output square matrix\

lns = cellfun(@length,map_in);
max_length = max(lns);
map_out = zeros(length(map_in),max_length);

for i = 1:length(map_in)
    map_out(i,1:lns(i)) = map_in{i};
end
