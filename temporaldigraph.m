function D = temporaldigraph(C,plotflag)

% Builds temporal digraph of a directed contact network in the sense of
% https://doi.org/10.1007/s41109-019-0209-1 (see here for definitions and
% code notation as well). No error checking is performed: caveat utilitor.
% 
% Example: C = [1,4,1;5,4,2;2,5,3;4,3,4]
%
% Last modified 20210906 by Steve Huntsman
%
% Copyright (c) 2021 Steve Huntsman
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

%% Eliminate any "loops" to make C a bona fide DCN
C = C(C(:,1)~=C(:,2),:);

%% Build fibers of DCN C
V = unique(C(:,[1,2]));
fiber = cell(1,numel(V));
for j = 1:numel(V)
    v = V(j);
    ind = union(find(C(:,1)==v),find(C(:,2)==v));
    % Usage of unique necessary here in case contacts with same source or
    % same target have identical timestamps
    fiber{j} = [-Inf;unique(C(ind,3));Inf];
end

%% Build temporal digraph D of DCN C
foo = [C(:,[1,3]),C(:,[2,3])];
bar = [];
for j = 1:numel(V)
    vn = repmat(V(j),size(fiber{j}));
    bar = [bar;[vn(1:end-1),fiber{j}(1:end-1),vn(2:end),fiber{j}(2:end)]];
end
AT = unique([foo;bar],'rows');
src = join(string(num2cell(AT(:,1:2))),",");
tar = join(string(num2cell(AT(:,3:4))),",");
D = digraph(src,tar);

%% Plot
if plotflag
    figure; 
    plo = plot(D);
    plo.NodeLabelMode = 'auto'; % force node labels
    nums = cell2mat(cellfun(@(s) sscanf(s,'%f,').',...
        plo.NodeLabel,'UniformOutput',false));
    nums = reshape(nums',[2,numel(nums)/2]);
    dt = min(diff(sort(C(:,3))));
    nums(nums==-Inf) = min(nums(2,isfinite(nums(2,:))))-dt;
    nums(nums==Inf) = max(nums(2,isfinite(nums(2,:))))+dt;
    plo.XData = nums(2,:);
    plo.YData = nums(1,:);
end