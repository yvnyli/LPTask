% variance (=???) vs confidence (=probability of decision being correct)

% confidence: e.g. given 3H1T, calculate probability of it coming from box
% of each mix of 100, normalize, the normalized prob from the mix of the
% decision (for example, the Bayesian decision or the MLE) is the
% confidence (meaning, it doesn't matter how you got to the decision, the
% confidence should be the same)


% variance: well i used the formula (meaning it only applies to Bayesian
% estimates)

% from the shape of the bottom two curves (confidence of bayesian estimates
% and inverse variance of bayesian posterior), i'm convinced that variance
% measures probability of being correct as well


n = 20;
cs = NaN(1,n+1);
for na = 0:n
  cs(na+1) = confidence(n,na,na/n,10000);
end
csbayes = NaN(1,n+1);
invvsbayes = NaN(1,n+1);
for na = 0:n
  csbayes(na+1) = confidence(n,na,(na+1)/(n+2),100);
  invvsbayes(na+1) = 1/(((na+1)*(n-na+1))/((n+2)^2*(n+3)));
end
figure;subplot(3,1,1);plot(cs);title('cs');
subplot(3,1,2);plot(csbayes);title('csbayes');
subplot(3,1,3);plot(invvsbayes);title('invvsbayes');
%%
function c = confidence(n,na,p,total)
if isempty(total)
  total = 100;
end
p = round(p * total)/total;
xs = NaN(1,total+1);
% work in log space, final answer = exp( logpprob - normalizer )
% where normalizer = x* + log(sum_i(exp(x_i)))
% where x_i = na * log(p_i) + (n-na) * log(1-p_i) and x* = max(xs)
if (p==0 && na==0) || (p==1 && na==n) % special case where we have 0*-Inf
  logpprob = 0;
else
  logpprob = na * log(p) + (n-na) * log(1-p);
end
for ind = 0:total
  p = ind/total;
  if (p==0 && na==0) || (p==1 && na==n) % special case where we have 0*-Inf
    xs(ind+1) = 0;
  else
    xs(ind+1) = na * log(p) + (n-na) * log(1-p);
  end
end
xmax = max(xs);
normalizer = xmax + log(sum(exp(xs-xmax)));
c = exp(logpprob-normalizer) * total;
end

% function c = confi(n,na,p,total)
% % i wonder if i need logsumexp at all
% if isempty(total)
%   total = 100;
% end
% p = round(p * total)/total;
% pprob = p^na * (1-p)^(n-na);
% normalizer = 0;
% for ind = 0:total
%   p = ind/total;
%   % if p=0 (p=1) then only if na=0 (na=n), prob is 1; otherwise prob is 0
%   normalizer = normalizer + p^na * (1-p)^(n-na);
% end
% c = pprob / normalizer;
% end
