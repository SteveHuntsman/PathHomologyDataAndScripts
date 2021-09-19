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
   C = [src(ind);tar(ind);tau(ind)-min(tau(ind))]';
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

%%
bettiEmail = betti;
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
filedir = '/Users/sha0639/Downloads/';
timestr = datestr(datetime,'yyyymmdd_HHMM');
filename = ['EmailBettiTemporalDigraph',timestr];
print('-dpdf',[filedir,filename,'.pdf'],'-r600');
print('-dpng',[filedir,filename,'.png'],'-r600');


