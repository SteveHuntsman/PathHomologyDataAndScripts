%% Get data from http://konect.uni-koblenz.de/networks/facebook-wosn-wall
fileDir = '/Users/sha0639/Downloads/';
filename = [fileDir,'out.facebook-wosn-wall.txt'];
str = fileread(filename);
str = strsplit(str);
raw = str2double(str(8:end));
src = raw(1:4:end);
tar = raw(2:4:end);
wgt = raw(3:4:end);
tau = raw(4:4:end);
N = min([numel(src),numel(tar),numel(tau)]);
src = src(1:N);
tar = tar(1:N);
wgt = wgt(1:N);
tau = tau(1:N);
clear raw str;

%%
[tau,ind] = sort(tau);
src = src(ind);
tar = tar(ind);
wgt = wgt(ind);
time = datetime(tau,'ConvertFrom','epochtime',...
'Format','dd-MMM-yyyy HH:mm:ss.SSS');

%% This has some obvious errors for j large; we recompute there symbolically
% The errors (in the form of negative Betti numbers) are caused by svds and
% can be fixed with symbolic calculations (albeit at large cost)
temp = datevec(time(1));
time0 = datetime(temp(1:3));
dt = 1;%8/24;
n_dt = 1600/dt; %
betti = zeros(n_dt,3);
%
for j = 1:1000%1:n_dt
   [j,n_dt]
   t1 = time0+j*dt+.3; % to capture a calendar day
   ind = isbetween(time,t1-dt,t1);
   % Avoid multiple edges and extraneous nodes
   ST = unique([src(ind);tar(ind)]','rows');
   D = digraph(ST(:,1),ST(:,2));
   D = subgraph(D,unique(D.Edges.EndNodes));
   D = removeloopsisos(D); % necessary here because of isolated loops
   %
   if size(D.Nodes,1)
       ph = pathhomology(D,size(betti,2));
       betti(j,:) = ph.betti;
   end
   %%
end

%% 
hasBneg = find(any(betti<0,2));

%% Recompute symbolically for anything with betti_2 purportedly negative
for i = 1:numel(hasBneg) 
   i
   j = hasBneg(i);
   t1 = time0+j*dt+.3; % to capture a calendar day
   ind = isbetween(time,t1-dt,t1);
   % Avoid multiple edges and extraneous nodes
   ST = unique([src(ind);tar(ind)]','rows');
   D = digraph(ST(:,1),ST(:,2));
   D = subgraph(D,unique(D.Edges.EndNodes));
   D = removeloopsisos(D); % necessary here because of isolated loops
   %
   if size(D.Nodes,1)
       ph = pathhomology(D,size(betti,2),'symbolic');
       betti(j,:) = ph.betti;
   end
end

%%
save('workspace20200717ph');

%%
j = find(betti(:,3)>0,1,'first');   % 756
t1 = time0+j*dt+.3; % to capture a calendar day
ind = isbetween(time,t1-dt,t1);
% Avoid multiple edges and extraneous nodes
ST = unique([src(ind);tar(ind)]','rows');
D = digraph(ST(:,1),ST(:,2));
D = subgraph(D,unique(D.Edges.EndNodes));
D = removeloopsisos(D); % necessary here because of isolated loops
% Visually inspect and isolate only connected component with beta_2 not
% obviously zero
figure; 
plot(D);
cc = conncomp(D,'Type','weak');
ind = find(cc==cc(192));
E = subgraph(D,ind);
%
figure;
plo = plot(E,'layout','force','EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
plo.NodeLabel = {};
ph = pathhomology(E,3,'representatives');
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
% filename = ['Facebook',timestr];
% print('-dpdf',[fileDir,filename,'.eps'],'-r600');
% print('-dpng',[fileDir,filename,'.png'],'-r600');
% Produce and save figure of all components 
figure; 
plo = plot(D,'EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
axis off;
hold on;
p = patch([0,5,5,0],[0,0,6,6],'k'); set(p,'FaceAlpha',0);
% filename = ['FacebookContext',timestr];
% print('-dpdf',[fileDir,filename,'.eps'],'-r600');
% print('-dpng',[fileDir,filename,'.png'],'-r600');

%% Plot Betti numbers and posts per day
h = figure; 
set(h,'Position',[0,0,300,200]);
temp = find(betti(1:1000,3));
hold on;
loglog(1:1000,betti(1:1000,1),'Marker','.','Color',0*[1,1,1],'LineStyle','none');
loglog(1:1000,betti(1:1000,2),'Marker','.','Color',.5*[1,1,1],'LineStyle','none');
temp = find(betti(1:1000,3));
loglog(temp,betti(temp,3),'Marker','o','Color',0*[1,1,1],'LineStyle','none');
set(gca,'XScale','log','YScale','log');
legend({'$\tilde \beta_0$','$\tilde \beta_1$','$\tilde \beta_2$'},...
   'Location','northwest','Interpreter','latex');
xlabel('days since beginning','Interpreter','latex');
box on;
ax = axis;
% filename = ['FacebookBetti',timestr];
% print('-dpdf',[fileDir,filename,'.eps'],'-r600');
% print('-dpng',[fileDir,filename,'.png'],'-r600');

%% Posts per day
ppd = zeros(1,1000);
for j = 1:1000
   t1 = time0+j*dt+.3; % to capture a calendar day
   ind = isbetween(time,t1-dt,t1);
   ppd(j) = nnz(ind);
end
h = figure;
set(h,'Position',[0,0,300,200]);
loglog(ppd,'k.');
xlabel('days since beginning','Interpreter','latex');
legend({'Posts per day'},...
   'Location','northwest','Interpreter','latex');
axis(ax);
box on;
% filename = ['FacebookActivity',timestr];
% print('-dpdf',[fileDir,filename,'.eps'],'-r600');
% print('-dpng',[fileDir,filename,'.png'],'-r600');
