%% NB. Various code cells need to be re-run with various filenames
% filename = '/Users/sha0639/Downloads/sx-mathoverflow-a2q.txt';      % http://snap.stanford.edu/data/sx-mathoverflow.html
% filename = '/Users/sha0639/Downloads/email-Eu-core-temporal.txt';   % http://snap.stanford.edu/data/email-Eu-core-temporal.html
% % % filename = '/Users/sha0639/Downloads/facebook-wall.txt';            % http://socialnetworks.mpi-sws.org/data-wosn2009.html      
filename = '/Users/sha0639/Downloads/out.facebook-wosn-wall.txt';   % Konect went down: used https://raw.githubusercontent.com/AbdelouahabKhelifati/EPredictor/88a08a9d4d109f4f73d5aaee35896b5032b000df/etc/facebook-wosn-wall/out.facebook-wosn-wall

%% Get temporal network in memory
if strcmp(filename,'/Users/sha0639/Downloads/out.facebook-wosn-wall.txt')
str = fileread(filename);
str = strsplit(str);
raw = str2double(str(8:end));
src = raw(1:4:end);
tar = raw(2:4:end);
tau = raw(4:4:end);
clear raw str;
else
A = readmatrix(filename);
src = A(:,1);
tar = A(:,2);
tau = A(:,3);
end
N = min([numel(src),numel(tar),numel(tau)]);
src = src(1:N);
tar = tar(1:N);
tau = tau(1:N);
[tau,ind] = sort(tau);
src = src(ind);
tar = tar(ind);
time = datetime(tau,'ConvertFrom','epochtime',...
'Format','dd-MMM-yyyy HH:mm:ss.SSS');

%% For node name prefixes (now doing this inline)
% Xsrc = cellfun(@num2str,num2cell(src),'UniformOutput',false);
% Xsrc = strcat(repmat({'x'},[numel(Xsrc),1]),Xsrc(:));
% Xtar = cellfun(@num2str,num2cell(tar),'UniformOutput',false);
% Xtar = strcat(repmat({'x'},[numel(Xtar),1]),Xtar(:));

%% Temporal parameters
if strcmp(filename,'/Users/sha0639/Downloads/sx-mathoverflow-a2q.txt')
    % After varying dt and window there is no evidence of homology in
    % dimension > 1 for any choice of these parameters, so using minimal
    % window in the interest of speed (by way of reproducibility)
    dt = 24/24;
    window = dt;
    T = days(time(end)-time(1));
    t0 = time(1);
elseif strcmp(filename,'/Users/sha0639/Downloads/email-Eu-core-temporal.txt')
    dt = 1/24;
    window = 2*dt;
    T = days(time(end)-time(1));
    t0 = time(1);
elseif strcmp(filename,'/Users/sha0639/Downloads/out.facebook-wosn-wall.txt')
    dt = 24/24;
    window = 2*dt;
    T = 1000;
    t0 = time(1)+.3;    % for boundaries in daily lulls
else
    error('what file?');
end
slide = dt;
t1 = t0+window;

%% Represent a la https://doi.org/10.1088/1367-2630/16/12/123055
dim = window/dt+1;
betti = [];
caught = [];
clustCoeff = [];
edgeDensity = [];
while t1<time(1)+T
    [t0,t1]
    D = digraph;
    for j = 1:window/dt
        ind01 = find(isbetween(time,t0+(j-1)*dt,t0+j*dt));
        for k = 1:numel(ind01)
            i = ind01(k);
            %             D = addedge(D,[Xsrc{i},'_',num2str(j)],[Xtar{i},'_',num2str(j+1)]);
            D = addedge(D,['x',num2str(src(i)),'_',num2str(j)],...
                ['x',num2str(tar(i)),'_',num2str(j+1)]);
        end
        D = simplify(D);
        %
        foo = localClustCoeff(adjacency(D));
        foo(isnan(foo)) = 0;
        clustCoeff = [clustCoeff;mean(foo)];
        edgeDensity = [edgeDensity;numedges(D)/(numnodes(D)*(numnodes(D)-1))];
    end
    if size(D.Nodes,1)
        try
            ph = pathhomology(D,dim);
            b01 = ph.betti;
        catch
            caught = [caught;t0];
            b01 = nan(1,dim);
        end
    else
        b01 = zeros(1,dim);
    end
    b01
    betti = [betti;b01];
    t0 = t0+slide;
    t1 = t0+window;
end
clustCoeff(isnan(clustCoeff)) = 0;
edgeDensity(isnan(edgeDensity)) = 0;

%% Vs
if strcmp(filename,'/Users/sha0639/Downloads/sx-mathoverflow-a2q.txt')
    figure('Position',[0,0,560,210]);
    subplot(1,2,1);
    plot(betti(:,2),clustCoeff,'k.');
    xlabel('$\tilde \beta_1$','Interpreter','latex');
    ylabel('clustering coefficient','Interpreter','latex');
    title('MathOverflow layered digraphs','Interpreter','latex');
    rho = corrcoef(betti(:,2),clustCoeff);
    legend(['$\rho = ',num2str(rho(1,2)),'$'],...
        'Interpreter','latex','Location','Northeast');
    subplot(1,2,2);
    plot(betti(:,2),edgeDensity,'k.');
    xlabel('$\tilde \beta_1$','Interpreter','latex');
    ylabel('edge density','Interpreter','latex');
    title('MathOverflow layered digraphs','Interpreter','latex');
    rho = corrcoef(betti(:,2),edgeDensity);
    legend(['$\rho = ',num2str(rho(1,2)),'$'],...
        'Interpreter','latex','Location','Northeast');
%     % Save figure to same directory
%     filedir = '/Users/sha0639/Downloads/';
%     timestr = datestr(datetime,'yyyymmdd_HHMM');
%     filename = ['VsMoMLP',timestr];
%     print('-dpdf',[fileDir,filename,'.pdf'],'-r600');
%     print('-dpng',[fileDir,filename,'.png'],'-r600');
elseif strcmp(filename,'/Users/sha0639/Downloads/email-Eu-core-temporal.txt')
    figure('Position',[0,0,560,210]);
    subplot(1,2,1);
    plot(betti(:,2),edgeDensity(1:2:end),'k.');
    xlabel('$\tilde \beta_1$','Interpreter','latex');
    ylabel('clustering coefficient','Interpreter','latex');
    title('European email layered digraphs','Interpreter','latex');
    rho = corrcoef(betti(:,2),edgeDensity(1:2:end));
    legend(['$\rho = ',num2str(rho(1,2)),'$'],...
        'Interpreter','latex','Location','Northeast');
    subplot(1,2,2);
    plot(betti(:,3),edgeDensity(1:2:end),'k.');
    xlabel('$\tilde \beta_2$','Interpreter','latex');
    ylabel('edge density','Interpreter','latex');
    title('European email layered digraphs','Interpreter','latex');
    rho = corrcoef(betti(:,3),edgeDensity(1:2:end));
    legend(['$\rho = ',num2str(rho(1,2)),'$'],...
        'Interpreter','latex','Location','Northeast');
%     % Save figure to same directory
%     filedir = '/Users/sha0639/Downloads/';
%     timestr = datestr(datetime,'yyyymmdd_HHMM');
%     filename = ['VsEmailMLP',timestr];
%     print('-dpdf',[fileDir,filename,'.pdf'],'-r600');
%     print('-dpng',[fileDir,filename,'.png'],'-r600');
elseif strcmp(filename,'/Users/sha0639/Downloads/out.facebook-wosn-wall.txt')
    figure('Position',[0,0,560,210]);
    subplot(1,2,1);
    plot(betti(:,2),edgeDensity(1:2:end),'k.');
    xlabel('$\tilde \beta_1$','Interpreter','latex');
    ylabel('clustering coefficient','Interpreter','latex');
    title('Facebook layered digraphs','Interpreter','latex');
    rho = corrcoef(betti(:,2),edgeDensity(1:2:end));
    legend(['$\rho = ',num2str(rho(1,2)),'$'],...
        'Interpreter','latex','Location','Northeast');
    subplot(1,2,2);
    plot(betti(:,3),edgeDensity(1:2:end),'k.');
    xlabel('$\tilde \beta_2$','Interpreter','latex');
    ylabel('edge density','Interpreter','latex');
    title('Facebook layered digraphs','Interpreter','latex');
    rho = corrcoef(betti(:,3),edgeDensity(1:2:end));
    legend(['$\rho = ',num2str(rho(1,2)),'$'],...
        'Interpreter','latex','Location','Northeast');
%     % Save figure to same directory
%     filedir = '/Users/sha0639/Downloads/';
%     timestr = datestr(datetime,'yyyymmdd_HHMM');
%     filename = ['VsEmailMLP',timestr];
%     print('-dpdf',[fileDir,filename,'.pdf'],'-r600');
%     print('-dpng',[fileDir,filename,'.png'],'-r600');
else
    error('what file?');
end

%% Just for saving stuff
if strcmp(filename,'/Users/sha0639/Downloads/sx-mathoverflow-a2q.txt')
bettiMO = betti;
elseif strcmp(filename,'/Users/sha0639/Downloads/email-Eu-core-temporal.txt')
bettiEmail = betti;
elseif strcmp(filename,'/Users/sha0639/Downloads/out.facebook-wosn-wall.txt')
bettiFB = betti;
else
error('what file?');
end

%% Plot Betti numbers for MO
h = figure('Position',[0,0,300,200]);
subplot(2,1,1);
loglog(1:size(bettiMO,1),bettiMO(:,1),'k.');
legend({'$\tilde \beta_0$'},...
'Location','northwest','Interpreter','latex');
axis([1,size(bettiMO,1),1,100]);
box on;
subplot(2,1,2);
loglog(1:size(bettiMO,1),bettiMO(:,2),'k.');
legend({'$\tilde \beta_1$'},...
'Location','northwest','Interpreter','latex');
axis([1,size(bettiMO,1),1,100]);
box on;
xlabel('days since beginning','Interpreter','latex');
% filedir = '/Users/sha0639/Downloads/';
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['MoBettiMLP',timestr];
% print('-dpdf',[filedir,filename,'.pdf'],'-r600');
% print('-dpng',[filedir,filename,'.png'],'-r600');

%% RERUN INITIAL MO STUFF BEFORE PLOTTING DIGRAPHS W/ 1-HOMOLOGY 
% There is no 2-homology (it turns out)
indMO2 = find(bettiMO(:,2)>20);
figure('Position',[0,0,560,420]);
ha = tight_subplot(3,1,.01,.01,.01);    % checked offline that this works
for im2 = 1:numel(indMO2)
t0 = time(1)+(indMO2(im2)-1)*slide; 
t1 = t0+window; 
D = digraph;
for j = 1:window/dt
   ind01 = find(isbetween(time,t0+(j-1)*dt,t0+j*dt));
   for k = 1:numel(ind01)
       i = ind01(k);
       D = addedge(D,['x',num2str(src(i)),'_',num2str(j)],...
           ['x',num2str(tar(i)),'_',num2str(j+1)]);
   end
   D = simplify(D);
end
% Decompose into weak components and analyze accordingly
cc = conncomp(D,'Type','weak');
ucc = unique(cc);
betticc = [];
phcc = cell(1,numel(ucc));
for j = 1:numel(ucc)
   G = subgraph(D,find(cc==ucc(j)));
   phcc{j} = pathhomology(G,dim,'representatives');
   betticc(j,:) = phcc{j}.betti;
end
keycc = find(betticc(:,end));
% Plot subgraph of each key weak component with highlighted representatives
   G = subgraph(D,ismember(cc,keycc));
   ph = pathhomology(G,dim,'representatives');
   [J,~] = find(ph.hom{dim});
   H = subgraph(G,unique(ph.allPaths{dim}(J,:)));
   ph = pathhomology(H,dim,'representatives'); % note reassignment
   [J,~] = find(ph.hom{dim});  % ibid.
   % h = figure;
   % set(h,'Position',[0,0,300,200]);
   axes(ha(im2));
   plo = plot(H,'layout','Layered','EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
%         for jj = 1:numel(J)
%            temp = ph.allPaths{dim}(J(jj),:);
%            sD = temp(1:end-1);
%            tD = temp(2:end);
%            highlight(plo,sD,tD,'EdgeColor',clr,'LineWidth',3);
%         end
   axis off;        
end
% filedir = '/Users/sha0639/Downloads/';
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['MoBetti1MLP',timestr];
% print('-dpdf',[filedir,filename,'.pdf'],'-r600');
% print('-dpng',[filedir,filename,'.png'],'-r600');    

%% Plot Betti numbers for enail
h = figure;
set(h,'Position',[0,0,560,200]);
loglog(1:size(bettiEmail,1),sort(bettiEmail(:,1),'descend'),'Marker','.','Color',0*[1,1,1],'LineStyle','none');
hold on;
loglog(1:size(bettiEmail,1),sort(bettiEmail(:,2),'descend'),'Marker','.','Color',.5*[1,1,1],'LineStyle','none');
temp = find(bettiEmail(:,3));
loglog(1:numel(temp),sort(bettiEmail(temp,3),'descend'),'Marker','o','Color',0*[1,1,1],'LineStyle','none');
set(gca,'XScale','log','YScale','log');
legend({'$\tilde \beta_0$','$\tilde \beta_1$','$\tilde \beta_2$'},...
'Location','Southwest','Interpreter','latex');
xlabel('rank','Interpreter','latex');
box on;
% filedir = '/Users/sha0639/Downloads/';
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['EmailBettiMLP',timestr];
% print('-dpdf',[filedir,filename,'.pdf'],'-r600');
% print('-dpng',[filedir,filename,'.png'],'-r600');

%% RERUN INITIAL EMAIL STUFF BEFORE PLOTTING DIGRAPHS W/ 2-HOMOLOGY 
% The array indEmail3 has repeated indices that are of interest as
% indicators of yet higher homology, but nothing pans out there. This plot
% is kind of a kludge. It shows, for each instance of Betti(2), the
% subgraph on vertices participating in 2-homology representatives. (It
% turns out in this particular instance that this subgraph is confined to a
% single weak component every time, which is fortunate or at least
% convenient since the code implicitly assumes that as written.)
indEmail3 = find(bettiEmail(:,3));
figure('Position',[0,0,560,700]);
ha = tight_subplot(6,3,.01,.01,.01);    % checked offline that this works
for ie3 = 1:numel(indEmail3)
t0 = time(1)+(indEmail3(ie3)-1)*slide; 
t1 = t0+window; 
D = digraph;
for j = 1:window/dt
   ind01 = find(isbetween(time,t0+(j-1)*dt,t0+j*dt));
   for k = 1:numel(ind01)
       i = ind01(k);
       D = addedge(D,['x',num2str(src(i)),'_',num2str(j)],...
           ['x',num2str(tar(i)),'_',num2str(j+1)]);
   end
   D = simplify(D);
end
% Decompose into weak components and analyze accordingly
cc = conncomp(D,'Type','weak');
ucc = unique(cc);
betticc = [];
phcc = cell(1,numel(ucc));
for j = 1:numel(ucc)
   G = subgraph(D,find(cc==ucc(j)));
   phcc{j} = pathhomology(G,dim,'representatives');
   betticc(j,:) = phcc{j}.betti;
end
keycc = find(betticc(:,end));
% Plot subgraph of each key weak component with highlighted representatives
for j = 1:numel(keycc)
   G = subgraph(D,find(cc==keycc(j)));
   ph = phcc{keycc(j)};
   [J,~] = find(ph.hom{dim});
   H = subgraph(G,unique(ph.allPaths{dim}(J,:)));
   ph = pathhomology(H,dim,'representatives'); % note reassignment
   [J,~] = find(ph.hom{dim});  % ibid.
   % h = figure;
   % set(h,'Position',[0,0,300,200]);
   axes(ha(ie3));
   plo = plot(H,'layout','Layered','EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
%         for jj = 1:numel(J)
%            temp = ph.allPaths{dim}(J(jj),:);
%            sD = temp(1:end-1);
%            tD = temp(2:end);
%            highlight(plo,sD,tD,'EdgeColor',clr,'LineWidth',3);
%         end
   axis off;        
end
end
% filedir = '/Users/sha0639/Downloads/';
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['EmailBetti2MLP',timestr];
% print('-dpdf',[filedir,filename,'.pdf'],'-r600');
% print('-dpng',[filedir,filename,'.png'],'-r600');

%% Plot Betti numbers for FB
h = figure;
set(h,'Position',[0,0,300,200]);
loglog(1:size(bettiFB,1),bettiFB(:,1),'Marker','.','Color',0*[1,1,1],'LineStyle','none');
temp = find(bettiFB(:,3));
hold on;
loglog(1:size(bettiFB,1),bettiFB(:,2),'Marker','.','Color',.5*[1,1,1],'LineStyle','none');
temp = find(bettiFB(:,3));
loglog(temp,bettiFB(temp,3),'Marker','o','Color',0*[1,1,1],'LineStyle','none');
set(gca,'XScale','log','YScale','log');
legend({'$\tilde \beta_0$','$\tilde \beta_1$','$\tilde \beta_2$'},...
'Location','northwest','Interpreter','latex');
xlabel('days since beginning','Interpreter','latex');
box on;
% filedir = '/Users/sha0639/Downloads/';
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['FacebookBettiMLP',timestr];
% print('-dpdf',[filedir,filename,'.pdf'],'-r600');
% print('-dpng',[filedir,filename,'.png'],'-r600');

%% RERUN INITIAL FB STUFF BEFORE PLOTTING DIGRAPH W/ 2-HOMOLOGY
% We need to have the FB temporal network in memory
indFB3 = find(bettiFB(:,3));
t0 = time(1)+.3+(indFB3(1)-1)*slide; 
t1 = t0+window; 
% Confirm Betti(2)
D = digraph;
for j = 1:window/dt
ind01 = find(isbetween(time,t0+(j-1)*dt,t0+j*dt));
for k = 1:numel(ind01)
   i = ind01(k);
%         D = addedge(D,[Xsrc{i},'_',num2str(j)],[Xtar{i},'_',num2str(j+1)]);
   D = addedge(D,['x',num2str(src(i)),'_',num2str(j)],...
       ['x',num2str(tar(i)),'_',num2str(j+1)]);
end
D = simplify(D);
end
ph = pathhomology(D,dim);
% Plotting D shows that node x1599_* is relevant: get its weak component.
% This could be automated by just working over all weak components
% separately, but we didn't do this.
cc = conncomp(D,'Type','weak');
key = cc(findnode(D,'x1599_1'));    % gives 41, as for _2, _3
G = subgraph(D,find(cc==key));
ph = pathhomology(G,dim,'representatives');
%
h = figure;
set(h,'Position',[0,0,300,200]);
plo = plot(G,'layout','Layered','EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
% plo.NodeLabel = {};
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
% filedir = '/Users/sha0639/Downloads/';
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['FacebookBetti2MLP',timestr];
% print('-dpdf',[filedir,filename,'.pdf'],'-r600');
% print('-dpng',[filedir,filename,'.png'],'-r600');

% %% Get data from http://konect.uni-koblenz.de/networks/facebook-wosn-wall
% filename = '/Users/sha0639/Downloads/out.facebook-wosn-wall.txt';
% str = fileread(filename);
% str = strsplit(str);
% raw = str2double(str(8:end));
% src = raw(1:4:end);
% tar = raw(2:4:end);
% wgt = raw(3:4:end);
% tau = raw(4:4:end);
% N = min([numel(src),numel(tar),numel(tau)]);
% src = src(1:N);
% tar = tar(1:N);
% wgt = wgt(1:N);
% tau = tau(1:N);
% clear raw str;
% 
% %%
% [tau,ind] = sort(tau);
% src = src(ind);
% tar = tar(ind);
% wgt = wgt(ind);
% time = datetime(tau,'ConvertFrom','epochtime',...
% 'Format','dd-MMM-yyyy HH:mm:ss.SSS');
% 
% %% This has some obvious errors for j large; we recompute there symbolically
% % The errors (in the form of negative Betti numbers) are caused by svds and
% % can be fixed with symbolic calculations (albeit at large cost)
% temp = datevec(time(1));
% time0 = datetime(temp(1:3));
% dt = 1;%8/24;
% n_dt = 1000/dt; %
% betti = zeros(n_dt,3);
% %
% for j = 1:1000%1:n_dt
%    [j,n_dt]
%    t1 = time0+j*dt+.3; % to capture a calendar day
%    ind = isbetween(time,t1-dt,t1);
%    % Avoid multiple edges and extraneous nodes
%    ST = unique([src(ind);tar(ind)]','rows');
%    D = digraph(ST(:,1),ST(:,2));
%    D = subgraph(D,unique(D.Edges.EndNodes));
%    D = removeloopsisos(D); % necessary here because of isolated loops
%    %
%    if size(D.Nodes,1)
%        ph = pathhomology(D,size(betti,2));
%        betti(j,:) = ph.betti;
%    end
%    %%
% end
% 
% %%
% hasBneg = find(any(betti<0,2));
% 
% %% Recompute symbolically for anything with betti_2 purportedly negative
% for i = 1:numel(hasBneg)
%    i
%    j = hasBneg(i);
%    t1 = time0+j*dt+.3; % to capture a calendar day
%    ind = isbetween(time,t1-dt,t1);
%    % Avoid multiple edges and extraneous nodes
%    ST = unique([src(ind);tar(ind)]','rows');
%    D = digraph(ST(:,1),ST(:,2));
%    D = subgraph(D,unique(D.Edges.EndNodes));
%    D = removeloopsisos(D); % necessary here because of isolated loops
%    %
%    if size(D.Nodes,1)
%        ph = pathhomology(D,size(betti,2),'symbolic');
%        betti(j,:) = ph.betti;
%    end
% end
% 
% %%
% save('workspace20200717ph');
% 
% %%
% j = find(betti(:,3)>0,1,'first');   % 756
% t1 = time0+j*dt+.3; % to capture a calendar day
% ind = isbetween(time,t1-dt,t1);
% % Avoid multiple edges and extraneous nodes
% ST = unique([src(ind);tar(ind)]','rows');
% D = digraph(ST(:,1),ST(:,2));
% D = subgraph(D,unique(D.Edges.EndNodes));
% D = removeloopsisos(D); % necessary here because of isolated loops
% % Visually inspect and isolate only connected component with beta_2 not
% % obviously zero
% figure;
% plot(D);
% cc = conncomp(D,'Type','weak');
% ind = find(cc==cc(192));
% E = subgraph(D,ind);
% %
% figure;
% plo = plot(E,'layout','force','EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
% plo.NodeLabel = {};
% [J,K] = find(ph.hom{3});
% for jj = 1:numel(J)
% %     if ph.hom{3}(J(jj),K(jj)) > 0, clr = [0,0,1]; else clr = [1,0,0]; end
%    clr = 0*[1,1,1];
%    temp = ph.allPaths{3}(J(jj),:);
%    sD = temp(1:end-1);
%    tD = temp(2:end);
%    highlight(plo,sD,tD,'EdgeColor',clr,'LineWidth',3);
% end
% axis off;
% % Save figure
% filedir = '/Users/sha0639/Downloads/';
% timestr = datestr(datetime,'yyyymmdd_HHMM');
% filename = ['Facebook',timestr];
% print('-dpdf',[filedir,filename,'.pdf'],'-r600');
% print('-dpng',[filedir,filename,'.png'],'-r600');
% % Produce and save figure of all components
% figure;
% plo = plot(D,'EdgeColor',0*[1,1,1],'NodeColor',0*[1,1,1]);
% axis off;
% hold on;
% p = patch([0,5,5,0],[0,0,6,6],'k'); set(p,'FaceAlpha',0);
% filename = ['FacebookContext',timestr];
% print('-dpdf',[filedir,filename,'.pdf'],'-r600');
% print('-dpng',[filedir,filename,'.png'],'-r600');
