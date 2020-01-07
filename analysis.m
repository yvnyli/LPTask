%%
% correlations: 
%  1st/2nd n and acc of 1st/2nd est,
%  1st/2nd n and var of 1st/2nd est, 
%  2nd-1st n and magnitude of delta est, 
%  delta est and H(true p),
%  magnitude of delta est and sign of delta H(est)
%% simulate n trials for a grid of p values with combos of 1st and 2nd obs
rng(0);
pgrid = 0.05:0.05:0.95;
nsims = 20;
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

%% box plot 1st/2nd n and acc of 1st/2nd est
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

%% box plot 1st/2nd n and var of 1st/2nd est
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

%% 2nd-1st n and magnitude of delta est
figure;
deltaest = squeeze(abs(trajectories(:,:,:,:,2)-trajectories(:,:,:,:,1)));
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