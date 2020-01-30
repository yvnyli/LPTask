% variance (=???) vs confidence (=probability of decision being correct)

% confidence: posterior probability (only applies to Bayesian estimates)


% variance: well i used the formula (only applies to Bayesian
% estimates)

% from the curves and table I would use confidence which is the same as
% mybetaposterior but in log space to compute the posterior probability
% which is by definition the probability of the decision being correct
% because even though betapdf gives the same (if divide by sampling rate
% i.e. 100) values, it is different close to 0 or 1


n = 20;
cs = NaN(1,n+1);
for na = 0:n
  cs(na+1) = betabinomialpmf(1+na,1+n-na,100,round(na/n*100));
end
csbayes = NaN(1,n+1);
csbayes2 = NaN(1,n+1);
csbayes3 = NaN(1,n+1);
invvsbayes = NaN(1,n+1);
for na = 0:n
  csbayes(na+1) = confidence(n,na,na/n,100);
  csbayes2(na+1) = betapdf(na/n,na+1,n-na+1);
  pmf = mybetaposterior(n,na,100);
  csbayes3(na+1) = pmf(round(na/n*100)+1);
  invvsbayes(na+1) = 1 / ( ((na+1)*(n-na+1)) / ((n+2)^2*(n+3)) );
end
figure;subplot(5,1,1);plot(cs);title('betabinomialpmf of mode');
xlim([1 21]);ylim([min(cs) max(cs)]);
subplot(5,1,2);plot(csbayes);title('confidence of mode');
xlim([1 21]);ylim([min(csbayes) max(csbayes)]);
subplot(5,1,3);plot(csbayes2);title('betapdf of mode');
xlim([1 21]);ylim([min(csbayes2) max(csbayes2)]);
subplot(5,1,4);plot(csbayes3);title('mybetaposterior of mode');
xlim([1 21]);ylim([min(csbayes3) max(csbayes3)]);
subplot(5,1,5);plot(invvsbayes);title('invvsbayes');
xlim([1 21]);ylim([min(invvsbayes) max(invvsbayes)]);
conftbl = table(cs',csbayes',csbayes2',csbayes3',invvsbayes')
%%
ns = [10,20,30,40,50,60,70,80,90];
cs = NaN(1,9);
for n = ns(:)'
  na = n/2;
  cs(n/10) = betabinomialpmf(1+na,1+n-na,100,round(na/n*100));
end
csbayes = NaN(1,9);
csbayes2 = NaN(1,9);
csbayes3 = NaN(1,9);
invvsbayes = NaN(1,9);
for n = ns(:)'
  na = n/2;
  csbayes(n/10) = confidence(n,na,na/n,100);
  csbayes2(n/10) = betapdf(na/n,na+1,n-na+1);
  pmf = mybetaposterior(n,na,100);
  csbayes3(n/10) = pmf(round(na/n*100)+1);
  invvsbayes(n/10) = 1 / ( ((na+1)*(n-na+1)) / ((n+2)^2*(n+3)) );
end
figure;subplot(5,1,1);plot(cs);title('betabinomialpmf of mode');
xlim([1 9]);ylim([min(cs) max(cs)]);
subplot(5,1,2);plot(csbayes);title('confidence of mode');
xlim([1 9]);ylim([min(csbayes) max(csbayes)]);
subplot(5,1,3);plot(csbayes2);title('betapdf of mode');
xlim([1 9]);ylim([min(csbayes2) max(csbayes2)]);
subplot(5,1,4);plot(csbayes3);title('mybetaposterior of mode');
xlim([1 9]);ylim([min(csbayes3) max(csbayes3)]);
subplot(5,1,5);plot(invvsbayes);title('invvsbayes');
xlim([1 9]);ylim([min(invvsbayes) max(invvsbayes)]);
conftbl2 = table(cs',csbayes',csbayes2',csbayes3',invvsbayes')
%%
function c = confidence(n,na,p,total)
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

function pmf = mybetaposterior(n,na,total)
pmf = NaN(1,total+1);
for ind = 0:total
  p = ind/total;
  pmf(ind+1) = p^na * (1-p)^(n-na);
end
pmf = pmf / sum(pmf);
end


function p = betabinomialpmf(a,b,n,k)
% the likelihood of observing k out of n samples given a beta distribution
p = nchoosek(n,k)*beta(k+a,n-k+b)/beta(a,b);
end

function c = confi(n,na,p,total)
% i wonder if i need logsumexp at all
if isempty(total)
  total = 100;
end
p = round(p * total)/total;
pprob = p^na * (1-p)^(n-na) * 1 / (total+1);
normalizer = 0;
for ind = 0:total
  p = ind/total;
  % if p=0 (p=1) then only if na=0 (na=n), prob is 1; otherwise prob is 0
  normalizer = normalizer + p^na * (1-p)^(n-na) * 1 / (total+1);
end
c = pprob / normalizer;
end
