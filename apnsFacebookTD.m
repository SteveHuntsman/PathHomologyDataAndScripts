%% Get data from http://konect.uni-koblenz.de/networks/facebook-wosn-wall
fileDir = '/Users/sha0639/Downloads/facebook-wosn-wall 2/';
filename = [fileDir,'out.facebook-wosn-wall'];
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

%% Path homology
% Any errors (in the form of negative Betti numbers) are caused by svds and
% can be fixed with symbolic calculations (albeit at large cost). We didn't
% see any here though.
temp = datevec(time(1));
time0 = datetime(temp(1:3));
dt = 1;%8/24;
n_dt = 1600/dt; %
betti = zeros(n_dt,2);
% Ran this with betti = zeros(n_dt,3) as well to see if beta_2 was ever
% nonzero: it eventually couldn't run because of too-large matrices, but
% beta_2 was never nonzero for the computed instances.
for j = 1:1000%1:n_dt
   [j,n_dt]
   t1 = time0+j*dt+.3; % to capture a calendar day
   ind = isbetween(time,t1-dt,t1);
   % Avoid multiple edges and extraneous nodes
   C = [src(ind);tar(ind);tau(ind)]';
   size(C,1)
   %
   if ~isempty(C)
       D = temporaldigraph(C,0);
       ph = pathhomology(D,size(betti,2),'removeLeaves')
       if any(ph.betti<0)
           ph = pathhomology(D,size(betti,2),'symbolic')
       end
       betti(j,:) = ph.betti;
   end
   %%
end

%% Plot Betti numbers and posts per day
h = figure; 
set(h,'Position',[0,0,300,200]);
hold on;
loglog(1:1000,betti(1:1000,1),'Marker','.','Color',0*[1,1,1],'LineStyle','none');
loglog(1:1000,betti(1:1000,2),'Marker','.','Color',.5*[1,1,1],'LineStyle','none');
set(gca,'XScale','log','YScale','log');
legend({'$\tilde \beta_0$','$\tilde \beta_1$'},...
   'Location','northwest','Interpreter','latex');
xlabel('days since beginning','Interpreter','latex');
box on;
ax = axis;
timestr = datestr(datetime,'yyyymmdd_HHMM');
filename = ['FacebookBettiTemporalDigraph',timestr];
print('-dpdf',[fileDir,filename,'.pdf'],'-r600');
print('-dpng',[fileDir,filename,'.png'],'-r600');

