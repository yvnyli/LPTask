pa = 0:0.1:1;
pb = 1-pa;
pc = max(pa,pb);
nsamples = 500;
updateweight = 100;
pcprime = pa.*( max( (pa*nsamples+updateweight)./(nsamples+updateweight),...
                     (pb*nsamples)./(nsamples+updateweight) )) + ...
  pb.*(max( (pa*nsamples)./(nsamples+updateweight),...
            (pb*nsamples+updateweight)./(nsamples+updateweight) ));
pg = pcprime-pc