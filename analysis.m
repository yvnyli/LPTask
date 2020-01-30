%%

% how to distangle the parameters ("model recovery"): simulate exp, suppose
% the agent generates the response (boredom) with a certain combination of
% the parameters (see google doc), fit the full regression model and see if
% the coefficients match the predetermined combination

% want a range for each of the parameters

%%%%%%%%%%%%%%%%%%%%%%%
%Boredom ~ b0, dConfidence, abs(dEstimate), conf1*abs(dEstimate), 
%          dEstimate, liking, n1, n2, est1, conf1

% surprise is a function of confidence1 and change in estimate
%   in the model, it can maybe be 1. the area of overlap between 1st and
%   2nd posterior dist., 2. log likelihood of 2nd subround of samples
%   based on 1st estimate, 3. 1/variance1 * absolute change in estimate


% bin abs(delta estimate) as conditions (say three) and have the same
% number of trials in each
% have all the combinations of n's ([4,8])
% and have these be kind of independent


%% all possible observations (n1,n2=4,8) and range for each of the parameters

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
% simulation
varNames = {'n1','n2','na1','na2','p1','adp','pp1','adpp',...
  'invvar1','dinvvar','negao','negll2','adppxinvvar1'};
tbl = table('Size',[(4+1)^2+(8+1)^2+(4+1)*(8+1)*2, numel(varNames)],...
  'VariableTypes',repmat({'double'},1,numel(varNames)),...
  'variablenames',varNames);
ns = [4,8];
trind = 0;
for n1 = ns(:)'
  for n2 = ns(:)'
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
%         tbl.e1(trind) = -p1*log2(p1)-(1-p1)*log2(1-p1);
%         if isnan(tbl.e1(trind)); tbl.e1(trind)=0; end
%         e2 = -p2*log2(p2)-(1-p2)*log2(1-p2);
%         if isnan(e2); e2=0; end
%         tbl.de(trind) = e2 - tbl.e1(trind);
        tbl.invvar1(trind) = 1/(((na1+1)*(n1-na1+1))/((n1+2)^2*(n1+3)));
        tbl.dinvvar(trind) = ...
          1/(((na1+na2+1)*(n1-na1+n2-na2+1))/((n1+n2+2)^2*(n1+n2+3))) - tbl.invvar1(trind);
        tbl.negao(trind) = 1-betaOverlap(na1+1,n1-na1+1,na1+na2+1,n1-na1+n2-na2+1);
        tbl.negll2(trind) = -(na2*log(pp1)+(n2-na2)*log(1-pp1)); 
        tbl.adppxinvvar1(trind) = tbl.adpp(trind)*tbl.invvar1(trind);
      end
    end
  end
end

%%%parameter correlations

ns = [4,8];
varnames = {'1st Est','abs d Est','1st inverse Var.','d inverse Var.','1-area','-logLikelihood','1/var x dEst'};
for n1 = ns(:)'
    for n2 = ns(:)'
        figure;
        for vind1 = 5:width(tbl)
            for vind2 = vind1:width(tbl)
                subplot(width(tbl)-4,width(tbl)-4,...
                    (width(tbl)-4)*(vind1-5)+vind2-4);
                scatter(tbl{tbl.n1==n1&tbl.n2==n2,vind2},...
                    tbl{tbl.n1==n1&tbl.n2==n2,vind1},10);
                ylabel(varnames{vind1-4});
                xlabel(varnames{vind2-4});
                sgtitle(sprintf('%d then %d',n1,n2))
            end
        end
    end
end

% surprise metrics only
ns = [4,8];
varnames = {'1-area','logLikelihood','1/var x dEst'};
for n1 = ns(:)'
    for n2 = ns(:)'
        figure;
        for vind1 = 1:3
            for vind2 = 1:3
                subplot(3,3,...
                    (3)*(vind1-1)+vind2);
                scatter(tbl{tbl.n1==n1&tbl.n2==n2,vind2+8},...
                    tbl{tbl.n1==n1&tbl.n2==n2,vind1+10},8);
                ylabel(varnames{vind1});
                xlabel(varnames{vind2});
                sgtitle(sprintf('%d then %d',n1,n2))
            end
        end
    end
end


%% %model recovery
%Boredom ~ b0, EST1, Surprize, Liking, n1, n2, trialNumber(timeontask) 
% EST1: p1 or pp1
% Surprise: adp, adpp, invvar1, dinvvar, negao, negll2, adppxinvvar1
% 14 combinations
models = cell(30,11,15);

% make a table of 60 trials
lmtbl = tbl(1:60,:); lmtbl.na1=[]; lmtbl.na2=[];
paramCombos = [1,1,0,0,0,0,0,0,0,0,0;...
               1,1,1,1,0,0,0,0,0,0,0;...
               1,1,1,0,0,1,0,0,0,0,0;...
               1,1,1,0,0,0,1,0,0,0,0;...
               1,1,1,0,0,0,0,1,0,0,0;...
               1,1,1,0,0,0,0,0,1,0,0;...
               1,1,1,0,0,0,0,0,0,1,0;...
               1,1,1,0,0,0,0,0,0,0,1;...
               1,1,0,1,1,0,0,0,0,0,0;...
               1,1,0,0,1,1,0,0,0,0,0;...
               1,1,0,0,1,0,1,0,0,0,0;...
               1,1,0,0,1,0,0,1,0,0,0;...
               1,1,0,0,1,0,0,0,1,0,0;...
               1,1,0,0,1,0,0,0,0,1,0;...
               1,1,0,0,1,0,0,0,0,0,1];
               
% test if a coefficient of 1 for each parameter can be recovered
for sind = 1:30 % subject ind
  for trind = 1:60 % trial ind
    lmtbl.n1(trind) = randi([6,10],1);
    lmtbl.n2(trind) = randi([6,10],1);
    na1 = randi([0,lmtbl.n1(trind)],1);
    lmtbl.p1(trind) = na1/lmtbl.n1(trind);
    na2 = randi([0,lmtbl.n2(trind)],1);
    lmtbl.adp(trind) = abs(na2/lmtbl.n2(trind) - lmtbl.p1(trind));
    lmtbl.pp1(trind) = (na1+1)/(lmtbl.n1(trind)+2);
    lmtbl.adpp(trind) = abs((1+na2)/(2+lmtbl.n2(trind)) - lmtbl.p1(trind));
    lmtbl.invvar1(trind) = 1/(((na1+1)*(lmtbl.n1(trind)-na1+1))/...
      ((lmtbl.n1(trind)+2)^2*(lmtbl.n1(trind)+3)));
    lmtbl.dinvvar(trind) = 1/(((na1+na2+1)*(lmtbl.n1(trind)-na1+lmtbl.n2(trind)-na2+1))/...
      ((lmtbl.n1(trind)+lmtbl.n2(trind)+2)^2*(lmtbl.n1(trind)+lmtbl.n2(trind)+3))) - lmtbl.invvar1(trind);
    lmtbl.negao(trind) = 1-betaOverlap(na1+1,lmtbl.n1(trind)-na1+1,na1+na2+1,...
      lmtbl.n1(trind)-na1+lmtbl.n2(trind)-na2+1);
    lmtbl.negll2(trind) = -(na2*log(pp1)+(lmtbl.n2(trind)-na2)*log(1-pp1)); 
    lmtbl.adppxinvvar1(trind) = lmtbl.adpp(trind)*lmtbl.invvar1(trind);
  end
  for pind = 1:11 % z-score each parameter over all trials for each subject
    lmtbl{:,pind} = (lmtbl{:,pind} - mean(lmtbl{:,pind})) / std(lmtbl{:,pind});
  end
  % fit one model as if the subject cares about 1 of each param at a time
  for pind = 1:11
    Y = 1 + 1 * lmtbl{:,pind};
    Y = Y + randn(size(Y))*0.5;
    lmtbl.Y = Y;
    for pcind = 1:15
      mdl = fitlm(lmtbl(:,logical([paramCombos(pcind,:),1])));
    %   plot(mdl)
      models{sind,pind,pcind} = mdl;
    end
  end
end
%%
clc
pind = 11;
% pcind = 1;
for pcind = 1:15
meanCoeff = 0;
isIncluded = paramCombos(pcind,pind);
sigCounter = table('size',[height(models{1,pind,pcind}.Coefficients),1],...
  'variabletypes',{'double'},'variablenames',{'sigCount'},'rownames',...
  models{1,pind,pcind}.Coefficients.Properties.RowNames);
sigCounter.sigCount = zeros([height(sigCounter),1]);
for sind = 1:30
%   disp(models{sind,pind,pcind})
  if isIncluded
    meanCoeff = meanCoeff + ...
      models{sind,pind,pcind}.Coefficients{lmtbl.Properties.VariableNames{pind},'Estimate'};
  end
  for tind = 1:height(sigCounter)
    sigCounter{tind,1} = sigCounter{tind,1} + double(...
      models{sind,pind,pcind}.Coefficients{tind,4}<0.05);
  end
end
meanCoeff = meanCoeff / 30; 
if isIncluded
  fprintf('\n\n\n\nMean coefficient of %s = %.3f\n\n',lmtbl.Properties.VariableNames{pind},meanCoeff);
else
  fprintf('\n\n\n\n %s is not in the model\n\n',lmtbl.Properties.VariableNames{pind});
end
disp(sigCounter);
end
% notes from 1/22
%  variance is not intuitive for confidence; confidence feels like n
%  boredom ~ 1, 1st Est, n1, n2, trialNumber(timeOnTask), surprise, liking
%  let's just randomly sample sequences, don't worry about cherry picking
%  trials
%  6-10 then 6-10
%  2 valences
%  parameter recovery: plot actual coeffs against recovered coeffs, measure
%  slop, r^2, etc
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
    histogram(tbl.dvar(tbl.n1==n1 & tbl.n2==n2),-0.04:0.0025:0.01,'orientation','horizontal'); 
    ylabel('number of sequences');xlabel('delta variance')
    title(sprintf('%d then %d',n1,n2));
  end
end
histedges = [0,0.07,0.7;...
  0,0.06,0.6;...
  -1,0.2,1;...
  -0.04,0.0025,0.01;...
  0,0.1,1;...
  -25,2.5,0];
for vind = 1:numel(rankvars)
  figure; % distribution of dvar: in most cases variance decreases (more confident)
  for ind1 = 1:3
    for ind2 = 1:3
      n1 = ns(ind1); n2 = ns(ind2); subplot(3,3,(ind1-1)*3+ind2);
      histogram(tbl{tbl.n1==n1 & tbl.n2==n2,rankvars{vind}},...
        linspace(histedges(vind,1),histedges(vind,3),15)); 
      ylabel('number of sequences');xlabel(rankvars{vind})
      title(sprintf('%d then %d',n1,n2));
    end
  end
end

% get rank or zscore (but they don't look normal) from distribution of parameters
rankvars = {'adp','adpp','de','dvar','ao','ll2'};
rankvarinds = [6,8,10,12,13,15];
tblrank = table('Size',[height(tbl), numel(rankvars)],...
  'VariableTypes',repmat({'double'},1,numel(rankvars)),...
  'variablenames',rankvars);
for n1 = [3,6,9]
  for n2 = [3,6,9]
    trinds = find(tbl.n1==n1&tbl.n2==n2);
    for vind = 1:numel(rankvars)
      c = unique(tbl{tbl.n1==n1&tbl.n2==n2,rankvars{vind}});
      for trind = trinds(:)'
        tblrank{trind,vind} = find(c==tbl{trind,rankvars{vind}},1);
      end
    end
  end
end
tblzscore = table('Size',[height(tbl), numel(rankvars)],...
  'VariableTypes',repmat({'double'},1,numel(rankvars)),...
  'variablenames',rankvars);
for n1 = [3,6,9]
  for n2 = [3,6,9]
    trinds = find(tbl.n1==n1&tbl.n2==n2);
    for vind = 1:numel(rankvars)
      tblzscore{trinds,vind} = zscore(tbl{trinds,rankvars{vind}},1);
    end
  end
end
%
figure; % compare distributions of adpp, dvar, and ll2
for trind = 1:height(tbl)
  ind1 = tbl.n1(trind)/3; ind2 = tbl.n2(trind)/3; subplot(3,3,(ind1-1)*3+ind2);
  hold on; plot([1 2 3],[tblrank.adpp(trind),tblrank.dvar(trind),tblrank.ll2(trind)],'o-');
end
for ind1 = 1:3
  for ind2 = 1:3
    n1 = ns(ind1); n2 = ns(ind2); subplot(3,3,(ind1-1)*3+ind2); 
    ylabel('rank (small to large)');xlabel('parameters');xticks([1 2 3]);xticklabels({'adpp','dvar','ll2'});
    title(sprintf('%d then %d',n1,n2));
  end
end
figure; % compare distributions of adpp, dvar, and ll2
for trind = 1:height(tbl)
  ind1 = tbl.n1(trind)/3; ind2 = tbl.n2(trind)/3; subplot(3,3,(ind1-1)*3+ind2);
  hold on; plot([1 2 3],[tblzscore.adpp(trind),tblzscore.dvar(trind),tblzscore.ll2(trind)],'o-');
end
for ind1 = 1:3
  for ind2 = 1:3
    n1 = ns(ind1); n2 = ns(ind2); subplot(3,3,(ind1-1)*3+ind2); 
    ylabel('z score');xlabel('parameters');xticks([1 2 3]);xticklabels({'adpp','dvar','ll2'});
    title(sprintf('%d then %d',n1,n2));
  end
end
%% surprise as area of overlap vs surprise as likelihood of 2nd given p1
a1 = 1; b1 = 10; a2 = 10; b2 = 1;
ao = betaOverlap(a1,b1,a2,b2)
%% test cases for surprise

% betaKL(11,11,1,1)/20 % should be less surprising than
% betaKL(11,11,6,6)/10 % failed
% 
% betaKL(11,11,6,6)/10 % should be less surprising than
% betaKL(16,6,6,6)/10 % success
% 
% betaKL(11,11,6,6)/10 % should be as surprising as or slightly more
% betaKL(21,21,6,6)/20 % failed
% 
% betaKL(6,16,1,11)/10 % should be more surprising than
% betaKL(1,21,1,11)/10 % success

% 1-betaOverlap(11,11,1,1) % should be less surprising than
% 1-betaOverlap(11,11,6,6) % failed
% 
% 1-betaOverlap(11,11,6,6) % should be less surprising than
% 1-betaOverlap(16,6,6,6) % success
% 
% 1-betaOverlap(11,11,6,6) % should be as surprising as or slightly more
% 1-betaOverlap(21,21,6,6) % failed
% 
% 1-betaOverlap(6,16,1,11) % should be more surprising than
% 1-betaOverlap(1,21,1,11) % success

% betaNegativeLogLikelihood(1,1,11,11) % should be less surprising than
% betaNegativeLogLikelihood(6,6,11,11) % failed
% 
% betaNegativeLogLikelihood(6,6,11,11) % should be less surprising than
% betaNegativeLogLikelihood(6,6,16,6) % success
% 
% betaNegativeLogLikelihood(6,6,11,11) % should be as surprising as or slightly more
% betaNegativeLogLikelihood(6,6,21,21) % failed
% 
% betaNegativeLogLikelihood(1,11,6,16) % should be more surprising than
% betaNegativeLogLikelihood(1,11,1,21) % success

abs(11/22-1/2)*betaInvVar(1,1) % should be less surprising than
abs(11/22-6/12)*betaInvVar(6,6) % failed: they are both 0

abs(11/22-6/12)*betaInvVar(6,6) % should be less surprising than
abs(16/22-6/12)*betaInvVar(6,6) % success

abs(11/22-6/12)*betaInvVar(6,6) % should be as surprising as or slightly more
abs(21/42-6/12)*betaInvVar(6,6) % success

abs(6/22-1/12)*betaInvVar(1,11) % should be more surprising than
abs(1/22-1/12)*betaInvVar(1,11) % success


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
function negll2 = betaNegativeLogLikelihood(a1,b1,a2,b2)
negll2 = -log(betabinomialpmf(a1,b1,a2+b2-a1-b1,a2-a1));
end
function p = betabinomialpmf(a,b,n,k)
% the likelihood of observing k out of n samples given a beta distribution
p = nchoosek(n,k)*beta(k+a,n-k+b)/beta(a,b);
end
function d = betaKL(pa,pb,qa,qb)
% p is the posterior, q is the prior, d measures the amount of information
% needed to update belief from q to p
klintegrand = @(x) betapdf(x,pa,pb) .* ...
  (log2(betapdf(x,pa,pb)) - log2(betapdf(x,qa,qb)));
d = integral(klintegrand,0,1);
end
function invvar = betaInvVar(a,b)
invvar = 1 / ((a*b)/((a+b)^2*(a+b+1)));
end
function c = betaProbability(n,na,p,total)
% calculate confidence as the posterior probability
% n: number of samples
% na: number of a's
% p: the probability will be of the estimate P(a) = p
% total: total number of cards in box (i.e. this is a discrete beta)
if isempty(total)
  total = 100;
end
p = round(p * total)/total;
xs = NaN(1,total+1);
% assume flat prior, meaning always use the total number of observations
% (1st + 2nd) to calculate confidence

% work in log space
logprior = log(1/(total+1));
if (p==0 && na==0) || (p==1 && na==n) % special case where we have 0*-Inf
  logpprob = log(1) + logprior;
else
  logpprob = na * log(p) + (n-na) * log(1-p) + logprior;
end
for ind = 0:total
  p = ind/total;
  if (p==0 && na==0) || (p==1 && na==n) % special case where we have 0*-Inf
    xs(ind+1) = log(1) + logprior;
  else
    xs(ind+1) = na * log(p) + (n-na) * log(1-p) + logprior;
  end
end
xmax = max(xs);
normalizer = xmax + log(sum(exp(xs-xmax)));
c = exp(logpprob-normalizer);
end
