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
betti = zeros(n_dt,2);
%
for j = 1:n_dt
   [j,n_dt]
   t1 = time0+j*dt;
   ind = isbetween(time,t1-3*dt,t1);
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

bettiMO = betti;

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
filedir = '/Users/sha0639/Downloads/';
timestr = datestr(datetime,'yyyymmdd_HHMM');
filename = ['MoBettiTemporalDigraph',timestr];
print('-dpdf',[filedir,filename,'.pdf'],'-r600');
print('-dpng',[filedir,filename,'.png'],'-r600');
