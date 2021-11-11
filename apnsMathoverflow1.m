%% Get data from http://snap.stanford.edu/data/sx-mathoverflow.html
fileDir = '/Users/sha0639/Downloads/';
filename = [fileDir,'sx-mathoverflow-a2q.txt'];
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
time = datetime(tau,'ConvertFrom','epochtime',...
   'Format','dd-MMM-yyyy HH:mm:ss.SSS');

%%
temp = datevec(time(1));
time0 = datetime(temp(1:3));
dt = 8/24;
n_dt = 2351/dt; % 2351 days
betti = zeros(n_dt,3);
clustCoeff = zeros(n_dt,1);
edgeDensity = zeros(n_dt,1);
%
for j = 1:n_dt
   % [j,n_dt]
   t1 = time0+j*dt;
   ind = isbetween(time,t1-3*dt,t1);
   % Avoid multiple edges and extraneous nodes
   ST = unique([src(ind);tar(ind)]','rows');
   D = digraph(ST(:,1),ST(:,2));
   D = subgraph(D,unique(D.Edges.EndNodes));
   %
   if size(D.Nodes,1)
       ph = pathhomology(D,size(betti,2));
       betti(j,:) = ph.betti;
       %
       foo = localClustCoeff(adjacency(D));
       foo(isnan(foo)) = 0;
       clustCoeff(j) = mean(foo);
       edgeDensity(j) = numedges(D)/(numnodes(D)*(numnodes(D)-1));
   end
   %%
end

%%
hasB2 = find(betti(:,3)>0);
j = hasB2(1);
t1 = time0+j*dt;
ind = find(isbetween(time,t1-3*dt,t1));
% Avoid multiple edges and extraneous nodes
[ST,indU] = unique([src(ind);tar(ind)]','rows');
assert(isequal(ST,[src(ind(indU))',tar(ind(indU))']),'ST error');
% Build digraph
temp = cellfun(@num2str,num2cell(unique(ST(:))),'UniformOutput',0);
NodeTable = cell2table(temp,'VariableNames',{'Name'});
temp = cellfun(@num2str,num2cell(ST),'UniformOutput',0);
EdgeTable = table(temp,time(ind(indU))','VariableNames',{'EndNodes','time'});
D = digraph(EdgeTable,NodeTable);
D = subgraph(D,unique(D.Edges.EndNodes));
% Compute path homology and plot
ph = pathhomology(D,size(betti,2),'representatives');
figure;
plo = plot(D,'layout','layered','EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
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
% filename = ['MathOverflow',timestr];
% print('-dpdf',[fileDir,filename,'.pdf'],'-r600');
% print('-dpng',[fileDir,filename,'.png'],'-r600');

%%
figure('Position',[0,0,560,210]);
plot(edgeDensity(hasB2),clustCoeff(hasB2),'ko',...
    edgeDensity,clustCoeff,'k.');
xlabel('clustering coefficient','Interpreter','latex');
ylabel('edge density','Interpreter','latex');
title('MathOverflow windowed digraphs','Interpreter','latex');
rho = corrcoef(betti(:,3),edgeDensity);
legend({'$\tilde \beta_2 > 0$'},...
    'Interpreter','latex','Location','Northeast');
% % Save figure to same directory
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['MathOverflowBettiVs',timestr];
% print('-dpdf',[fileDir,filename,'.pdf'],'-r600');
% print('-dpng',[fileDir,filename,'.png'],'-r600');