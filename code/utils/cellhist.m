function [n cellout]=cellhist(CELL)
% this function plots a histogram based on char cell array.
% Input: CELL - a cell string array (Nx1)
% 
%Output: n - alements in a bin
%        cellout - the bin value( a char)

%Example

if size(CELL,2)>1
    error('CELL need to be a vector of Nx1')
end

if sum(cellfun(@ischar,CELL))~=size(CELL,1)
    error('CELL must be a cell string array') 
end

[cellout, mm, nn] = unique(CELL);
for i=1:length(cellout)
    n(i,1)=sum(nn==i);
end
[n,IX] = sort(n);
cellout=cellout(IX);
bar(1:length(n),n);
1;
set(gca,'XTick',1:length(n))
set(gca,'XTickLabel',cellout)
rotateticklabel(gca,90)
