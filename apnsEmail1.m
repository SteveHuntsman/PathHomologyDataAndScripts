%% Get data from http://snap.stanford.edu/data/email-Eu-core-temporal.html
fileDir = '/Users/sha0639/Downloads/';
filename = [fileDir,'email-Eu-core-temporal.txt'];
str = fileread(filename);
raw = str2double(strsplit(str));
src = raw(1:3:end-1);
tar = raw(2:3:end);
tau = raw(3:3:end);
N = min([numel(src),numel(tar),numel(tau)]);
src = src(1:N);
tar = tar(1:N);
tau = tau(1:N);
clear raw str;

%%
[tau,ind] = sort(tau);
src = src(ind);
tar = tar(ind);

%%
n = 100;
dn = 50;
n_dn = ceil(N/dn);
betti = zeros(n_dn,3);
%
for j = 1:n_dn
   [j,n_dn]
   ind = ((j-1)*dn)+(1:n);
   % Avoid multiple edges and extraneous nodes
    try
       ST = unique([src(ind);tar(ind)]','rows')+1; % +1!
       D = digraph(ST(:,1),ST(:,2));
       D = subgraph(D,unique(D.Edges.EndNodes));
       %
       if size(D.Nodes,1)
           ph = pathhomology(D,size(betti,2));
           betti(j,:) = ph.betti;
       end
    catch
   %%
end

%%
hasB2 = find(betti(:,3)>20);
j = hasB2(2);
ind = ((j-1)*dn)+(1:n);
% Avoid multiple edges and extraneous nodes
[ST,indU] = unique([src(ind);tar(ind)]','rows');
assert(isequal(ST,[src(ind(indU))',tar(ind(indU))']),'ST error');
% Build digraph
temp = cellfun(@num2str,num2cell(unique(ST(:))),'UniformOutput',0);
NodeTable = cell2table(temp,'VariableNames',{'Name'});
temp = cellfun(@num2str,num2cell(ST),'UniformOutput',0);
EdgeTable = table(temp,tau(ind(indU))','VariableNames',{'EndNodes','time'});
D = digraph(EdgeTable,NodeTable);
D = subgraph(D,unique(D.Edges.EndNodes));
% Compute path homology and plot
ph = pathhomology(D,size(betti,2),'representatives');
figure;
plo = plot(D,'layout','layered','EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
plo.NodeLabel = {};
daspect([1,.1,1]);
[J,K] = find(ph.hom{3});
for jj = 1:numel(J)
%     if ph.hom{3}(J(jj),K(jj)) > 0, clr = [0,0,1]; else clr = [1,0,0]; end
   clr = 0*[1,1,1];
   temp = ph.allPaths{3}(J(jj),:);
   sD = temp(1:end-1);
   tD = temp(2:end);
   highlight(plo,sD,tD,'EdgeColor',clr,'LineWidth',3);
end
axis off;
% % Save figure to same directory
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['EuropeanEmail',timestr];
% print('-dpdf',[fileDir,filename,'.eps'],'-r600');
% print('-dpng',[fileDir,filename,'.png'],'-r600');

h = figure; 
set(h,'Position',[0,0,560,210]);
histogram(betti(:,3),'BinMethod','integers','FaceColor',.75*[1,1,1]);
set(gca,'yscale','log');
ylim([.5,1e4]);
xlabel('$\tilde \beta_2$','Interpreter','latex');
title('Frequency of $\tilde \beta_2$ in European email windowed digraphs','Interpreter','latex');
% % Save figure to same directory
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['EuropeanEmailHistogramBetti2',timestr];
% print('-dpdf',[fileDir,filename,'.eps'],'-r600');
% print('-dpng',[fileDir,filename,'.png'],'-r600');
