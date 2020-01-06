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

%%
colors = viridis(18);
figure;
for ind1 = 1:numel(nsamples)
  for ind2 = 1:numel(nsamples)
    for inds = 1:nsims
      plot([nsamples(ind1),nsamples(ind1)+nsamples(ind2)],...
        squeeze(trajectories(indp,inds,ind1,ind2,:))-trajectories(indp,inds,ind1,ind2,1),'o-',...
        'color',colors(9+nsamples(ind2)-nsamples(ind1),:));
      hold on
    end
  end
end
line([3,18],[0,0],'color','k');