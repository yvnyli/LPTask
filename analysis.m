%%

% how to distangle the parameters ("model recovery"): simulate exp, suppose
% the agent generates the response (boredom) with a certain combination of
% the parameters (see google doc), fit the full regression model and see if
% the coefficients match the predetermined combination

% want a range for each of the parameters

%%%%%%%%%%%%%%%%%%%%%%%
%Boredom ~ b0, dConfidence, abs(dEstimate), conf1*abs(dEstimate), 
%          dEstimate, liking, n1, n2, est1, conf1, entropy1, dentropy

% surprise is a function of confidence1 and change in estimate
%   in the model, it can maybe be the overlap of the distribution of the
%   estimate based on the first samples and the distribution of the estimate
%   based on the second samples (think hypothesis testing, confidence interval etc)

% OK, both posteriors are beta distributions, let's try to use the area of
% overlap as a metric. max(AO) = 1 when the two are identical, min(AO) = 0
% when the two are completely dirac deltas (doesn't matter where though)


% bin abs(delta estimate) as conditions (say three) and have the same
% number of trials in each
% have all the combinations of n's ([3 6 9])
% and have these be kind of independent

% play with the true p but hey its not a regressor and nobody cares
% but i still want to do the analysis where i generate a bunch of sampling
% sequences and select them according to the conditions and see whether
% they look more like sequences from a distribution with a different p
% or that they look more like sequences from a volatile env


%% all possible observations (n1,n2=3,6,9) and range for each of the parameters

% 3 then 3: 2^6 sequences, 4*4=16 count combinations
% 3 then 6 or 6 then 3: 2^9 sequences, 4*7=28 count combinations
% 3 then 9 or 9 then 3: 2^12, 4*10=40 count combinations
% 6 then 6: 2^12, 7*7=49 count combinations
% 6 then 9 or 9 then 6: 2^15, 7*10=70 count combinations
% 9 then 9: 2^18, 10*10=100 count combinations

% participant parameters:
%  initial confidence and delta confidence
%  surprise (confidence * absolute change in estimate)
%  liking

% model parameters:
%  initial posterior variance and delta posterior variance
%  distribution overlap or likelihood of 2nd subround given estimate 1

% shared parameters:
%  n1, change in n2
%  estimate 1, absolute change in estimate
%  entropy of estimate 1, change in entropy 

% p and e are frequentist, pp and ll2 are with prior of 1 and 1
% var and ao are also of course using Bayesian estimates
varNames = {'n1','n2','na1','na2','p1','adp','pp1','adpp','e1','de',...
  'var1','dvar','ao','ll2'};
tbl = table('Size',[16+28*2+40*2+49+70*2+100, numel(varNames)],...
  'VariableTypes',repmat({'double'},1,numel(varNames)),...
  'variablenames',varNames);
trind = 0;
for n1 = [3,6,9]
  for n2 = [3,6,9]
    for na1 = 0:n1
      for na2 = 0:n2
        trind = trind + 1;
        tbl.n1(trind) = n1; tbl.n2(trind) = n2; 
        tbl.na1(trind) = na1; tbl.na2(trind) = na2;
        p1 = na1/n1; tbl.p1(trind) = p1; 
        p2 = (na1+na2)/(n1+n2); tbl.adp(trind) = abs(p2-p1);
        pp1 = (na1+1)/(n1+2);
        tbl.pp1(trind) = pp1; 
        tbl.adpp(trind) = abs((na1+na2+1)/(n1+n2+2)-pp1);
        tbl.e1(trind) = -p1*log2(p1)-(1-p1)*log2(1-p1);
        if isnan(tbl.e1(trind)); tbl.e1(trind)=0; end
        e2 = -p2*log2(p2)-(1-p2)*log2(1-p2);
        if isnan(e2); e2=0; end
        tbl.de(trind) = e2 - tbl.e1(trind);
        tbl.var1(trind) = ((na1+1)*(n1-na1+1))/((n1+2)^2*(n1+3));
        tbl.dvar(trind) = ...
          ((na1+na2+1)*(n1-na1+n2-na2+1))/((n1+n2+2)^2*(n1+n2+3)) - tbl.var1(trind);
        tbl.ao(trind) = betaOverlap(na1+1,n1-na1+1,na1+na2+1,n1-na1+n2-na2+1);
        tbl.ll2(trind) = na2*log(pp1)+(n2-na2)*log(1-pp1); 
      end
    end
  end
end


%% range for each of the parameters
figure;scatter(tbl.ao,tbl.ll2,[],tbl.de);xlabel('area overlap');ylabel('log likelihood 2');title('colored by delta entropy');colorbar;colormap(viridis(100));
figure;scatter(tbl.ao,tbl.ll2,[],tbl.dvar);xlabel('area overlap');ylabel('log likelihood 2');title('colored by delta variance');colorbar;colormap(viridis(100));
figure;scatter(tbl.de,tbl.dvar,[],tbl.n2-tbl.n1);xlabel('delta entropy');ylabel('delta variance');title('colored by n2-n1');colorbar;colormap(viridis(100));

color369 = zeros(9,3);color369(9,:) = [0.9,0.05,0.05];color369(6,:) = [0.05,0.9,0.05];color369(3,:) = [0.05,0.05,0.9];
figure;
for n1 = [3,6,9]
  for n2 = [3,6,9]
    subplot(2,2,1);scatter(tbl.ao(tbl.n1==n1 & tbl.n2==n2),tbl.dvar(tbl.n1==n1 & tbl.n2==n2),20,'linewidth',1.5,'markeredgecolor',color369(n2,:),'markerfacecolor',color369(n1,:));xlabel('area overlap');ylabel('delta variance');title('ao vs dvar');hold on;
    subplot(2,2,2);scatter(tbl.ll2(tbl.n1==n1 & tbl.n2==n2),tbl.dvar(tbl.n1==n1 & tbl.n2==n2),20,'linewidth',1.5,'markeredgecolor',color369(n2,:),'markerfacecolor',color369(n1,:));xlabel('log likelihood 2');ylabel('delta variance');title('ll2 vs dvar');hold on;
    subplot(2,2,3);scatter(tbl.ao(tbl.n1==n1 & tbl.n2==n2),tbl.de(tbl.n1==n1 & tbl.n2==n2),20,'linewidth',1.5,'markeredgecolor',color369(n2,:),'markerfacecolor',color369(n1,:));xlabel('area overlap');ylabel('delta entropy');title('ao vs de');hold on;
    subplot(2,2,4);scatter(tbl.ll2(tbl.n1==n1 & tbl.n2==n2),tbl.de(tbl.n1==n1 & tbl.n2==n2),20,'linewidth',1.5,'markeredgecolor',color369(n2,:),'markerfacecolor',color369(n1,:));xlabel('log likelihood 2');ylabel('delta entropy');title('ll2 vs de');hold on;
  end
end

ns = [3,6,9];

figure; % distribution of dvar: in most cases variance decreases (more confident)
for ind1 = 1:3
  for ind2 = 1:3
    n1 = ns(ind1); n2 = ns(ind2); subplot(3,3,(ind1-1)*3+ind2);
    histogram(tbl.dvar(tbl.n1==n1 & tbl.n2==n2),-0.04:0.0025:0.01); 
    ylabel('number of sequences');xlabel('delta variance')
    title(sprintf('%d then %d',n1,n2));
  end
end
%% surprise as area of overlap vs surprise as likelihood of 2nd given p1
a1 = 1; b1 = 10; a2 = 10; b2 = 1;
ao = betaOverlap(a1,b1,a2,b2)
%%
% correlations: 
%  1st/2nd n and acc of 1st/2nd est,
%  1st/2nd n and var of 1st/2nd est, 
%  2nd-1st n and magnitude of delta est, 
%  magnitude of delta est and H(true p),
%  magnitude of delta est and sign of delta H(est)

% simulate n trials for a grid of p values with combos of 1st and 2nd obs
rng(0);
pgrid = 0.05:0.05:0.95;
nsims = 150;
nsamples = [3,5,7,9];
observations = NaN(numel(pgrid),nsims,numel(nsamples),numel(nsamples),2*max(nsamples));
trajectories = NaN(numel(pgrid),nsims,numel(nsamples),numel(nsamples),2);
for indp = randperm(numel(pgrid))
  p = pgrid(indp);
  for inds = 1:nsims
    for ind1 = 1:numel(nsamples)
      for ind2 = 1:numel(nsamples)
       obs1 = round(rand(nsamples(ind1),1)+p-0.5);
       observations(indp,inds,ind1,ind2,1:numel(obs1)) = obs1;
       trajectories(indp,inds,ind1,ind2,1) = ...
         sum(obs1)/numel(obs1);
       obs2 = round(rand(nsamples(ind2),1)+p-0.5);
       observations(indp,inds,ind1,ind2,numel(obs1)+(1:numel(obs2))) = obs2;
       trajectories(indp,inds,ind1,ind2,2) = ...
         (sum(obs1)+sum(obs2))/(numel(obs1)+numel(obs2));
      end
    end
  end
end


% box plot 1st/2nd n and acc of 1st/2nd est
% acc = true p - est
% in exp, acc = participant est - bayesian est

figure;
% 1st
t1st = squeeze(trajectories(:,:,:,:,1)) - repmat(pgrid',1,nsims,numel(nsamples),numel(nsamples));
n1st = NaN(size(t1st));
for indp = 1:numel(pgrid)
  for inds = 1:nsims
    n1st(indp,inds,:,:) = repmat(nsamples',1,numel(nsamples));
  end
end
subplot(2,2,1);
boxplot(t1st(:),n1st(:))
set(gca,'xticklabels',num2str(nsamples'));
title('true p - bayesian est by number of obs (first)');

% 2nd
t2st = squeeze(trajectories(:,:,:,:,2)) - repmat(pgrid',1,nsims,numel(nsamples),numel(nsamples));
n2st = NaN(size(t2st));
for indp = 1:numel(pgrid)
  for inds = 1:nsims
    for ind1 = 1:numel(nsamples)
      n2st(indp,inds,ind1,:) = nsamples;
    end
  end
end
subplot(2,2,2);
boxplot(t2st(:),n2st(:))
set(gca,'xticklabels',num2str(nsamples'));
title('true p - bayesian est by number of obs (second)');


% box plot 1st/2nd n and var of 1st/2nd est
% var over sims only

% 1st
var1 = squeeze(std(trajectories(:,:,:,:,1),0,2));
n1st = NaN(size(var1));
for indp = 1:numel(pgrid)
  n1st(indp,:,:) = repmat(nsamples',1,numel(nsamples));
end
subplot(2,2,3);
boxplot(var1(:),n1st(:))
set(gca,'xticklabels',num2str(nsamples'));
title('variance of bayesian est by number of obs (first)');

% 2nd
var2 = squeeze(std(trajectories(:,:,:,:,2),0,2));
n2st = NaN(size(var2));
for indp = 1:numel(pgrid)
  for ind1 = 1:numel(nsamples)
    n2st(indp,ind1,:) = nsamples;
  end
end
subplot(2,2,4);
boxplot(var2(:),n2st(:))
set(gca,'xticklabels',num2str(nsamples'));
title('variance of bayesian est by number of obs (second)');


% 2nd-1st n and magnitude of delta est
figure;
% should i sqrt to make it more gaussian?
%  wikipedia says "If the sample size is large and the population is not 
%  normal, then the sample correlation coefficient remains approximately 
%  unbiased, but may not be efficient.
deltaest = squeeze((abs(trajectories(:,:,:,:,2)-trajectories(:,:,:,:,1))));
n2m1 = NaN(size(deltaest));
for indp = 1:numel(pgrid)
  for inds = 1:nsims
    for ind1 = 1:numel(nsamples)
      for ind2 = 1:numel(nsamples)
        n2m1(indp,inds,ind1,ind2) = ind2-ind1;
      end
    end
  end
end
boxplot(deltaest(:),n2m1(:));
title('abs change in estimate by 2nd n - 1st n')
% figure;
% subplot(3,1,1);normplot(deltaest(:));title('deltaest');
% subplot(3,1,2);normplot(sqrt(deltaest(:)));title('sqrt');
% subplot(3,1,3);normplot(log(deltaest(:)+0.01));title('log+0.01');
% deltaest = deltaest(:); n2m1 = n2m1(:);
% [r,p] = corr(deltaest(deltaest>0),n2m1(deltaest>0))

% magnitude of delta est and H(true p)
figure;
deltaest = squeeze(abs(trajectories(:,:,:,:,2)-trajectories(:,:,:,:,1)));
hp = NaN(size(deltaest));
truep = NaN(size(hp));
for indp = 1:numel(pgrid)
  p = pgrid(indp);
  hp(indp,:,:,:) = repmat(-log(p)*p-log(1-p)*(1-p),1,size(hp,2),size(hp,3),size(hp,4));
  truep(indp,:,:,:) = repmat(p,1,size(truep,2),size(truep,3),size(truep,4));
end
mdl = fitlm(hp(:),deltaest(:),'varnames',{'entropy_of_true_p','abs_change_in_estimate'})
plot(mdl);

figure;
for ind1 = 1:4
  for ind2 = 1:4
    subplot(4,4,ind2+(ind1-1)*4);
    deltaest = squeeze(trajectories(:,:,ind1,ind2,2)-trajectories(:,:,ind1,ind2,1));
    hp = NaN(size(deltaest));
    for indp = 1:numel(pgrid)
      p = pgrid(indp);
      hp(indp,:) = repmat(-log(p)*p-log(1-p)*(1-p),1,size(hp,2));
%       deltaest(indp,:) = (deltaest(indp,:) - mean(deltaest(indp,:),'all')) ./ std(deltaest(indp,:),0,'all');
    end
    mdl = fitlm(hp(:),abs(deltaest(:)),'varnames',{'entropy_of_true_p','abs_change_in_estimate'});
    plot(mdl);legend('off');
    title([num2str(nsamples(ind1)) ' then ' num2str(nsamples(ind2)) ', b=',...
      num2str(mdl.Coefficients.Estimate(2)),', p=',num2str(mdl.Coefficients.pValue(2))]);
  end
end


% magnitude of delta est and sign of delta H(est)
figure;
deltaest = squeeze((trajectories(:,:,:,:,2)-trajectories(:,:,:,:,1)));
deltah = NaN(size(deltaest));
for indp = 1:numel(pgrid)
  for inds = 1:nsims
    for ind1 = 1:numel(nsamples)
      for ind2 = 1:numel(nsamples)
        e1 = trajectories(indp,inds,ind1,ind2,1);
        e1 = max(min(e1,0.99),0.01);
        e2 = trajectories(indp,inds,ind1,ind2,2);
        e1 = max(min(e1,0.99),0.01);
        deltah(indp,inds,ind1,ind2) = (-log(e2)*e2-log(1-e2)*(1-e2)) - ...
          (-log(e1)*e1-log(1-e1)*(1-e1));
      end
    end
  end
end
mdl = fitlm(deltah(:),deltaest(:),'varnames',{'change_in_entropy','abs_change_in_estimate'});
plot(mdl);legend('off');















%% functions
function ao = betaOverlap(a1,b1,a2,b2)
overlap = @(x) min(betapdf(x,a1,b1),betapdf(x,a2,b2));
ao = integral(overlap,0,1);
end
