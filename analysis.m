%%
% correlations: 
%  1st/2nd n and acc of 1st/2nd est,
%  1st/2nd n and var of 1st/2nd est, 
%  2nd-1st n and magnitude of delta est, 
%  magnitude of delta est and H(true p),
%  magnitude of delta est and sign of delta H(est)
%% simulate n trials for a grid of p values with combos of 1st and 2nd obs
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
%% magnitude of delta est and H(true p)
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


%% magnitude of delta est and sign of delta H(est)
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