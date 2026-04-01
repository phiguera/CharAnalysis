function [mu,sig,ll,prop_n,optmixture] = GaussianMixture(pixels, initK,...
    finalK,verbose)
% [mixture, optmixture] = GaussianMixture(pixels, initK, finalK)
%        perform the EM algorithm to estimate the order, and
%        parameters of a Gaussian Mixture model for a given set of
%        observation.
%
%     pixels: a N x M matrix of observation vectors with each row being an
%        M-dimensional observation vector, totally N observations
%     initK: the initial number of clusters to start with and will be reduced
%        to find the optimal order or the desired order based on MDL
%     finalK: the desired final number of clusters for the model.
%        Estimate the optimal order if finalK == 0.
%     verbose: true/false, return clustering information if true
%
%     mixture: an array of mixture structures with each containing the
%        converged Gaussian mixture at a given order
%           mixture(l).K: order of the mixture
%           mixture(l).M: dimension of observation vectors
%           mixture(l).rissanen: converaged MDL(K)
%           mixture(l).loglikelihood: ln( Prob{Y=y|K, theta*} )
%           mixture(l).cluster: an array of cluster structures with each
%              containing the converged cluster parameters
%                 mixture(l).cluster(k).pb: pi(k)=Prob(Xn=k|K, theta*)
%                 mixture(l).cluster(k).mu: mu(k), mean vector of the k-th cluster
%                 mixture(l).cluster(k).R: R(k), covariance matrix of the k-th cluser
%     optmixture: one of the element in the mixture array.
%        If finalK > 0, optmixture = mixture(1) and is the mixture with order finalK.
%        If finalK == 0, optmixture is the one in mixture with minimum MDL
%
%   v2.0: Supporting functions (EStep, MStep, EMIterate, MDLReduceOrder,
%   ClusterNormalize, initMixture, SplitClasses, GMClassLikelihood) have
%   been moved from separate .m files into local functions within this file.
%   All code is verbatim from the original Bowman CLUSTER implementation
%   as distributed with CharAnalysis v1.1 — no logic changes.

if ~isnumeric(initK) || ~all(size(initK)==[1,1]) || initK<=0 || mod(initK,1)~=0
   error('GaussianMixture: initial number of clusters initK must be a positive integer');
end
if ~isnumeric(finalK) || ~all(size(finalK)==[1,1]) || finalK<0 || mod(finalK,1)~=0
   error('GaussianMixture: final number of clusters finalK must be a positive integer or zero');
end
if finalK > initK
   error('GaussianMixture: finalK cannot be greater than initK');
end
if ~isa(pixels,'double')
   pixels = double(pixels);
end
mtr = initMixture(pixels, initK);
mtr = EMIterate(mtr, pixels);
if verbose
   disp(sprintf('K: %d\t rissanen: %f', mtr.K, mtr.rissanen));
end
mixture(mtr.K-max(1,finalK)+1) = mtr;
while mtr.K > max(1, finalK)
   mtr = MDLReduceOrder(mtr, verbose);
   mtr = EMIterate(mtr, pixels);
   if verbose
      disp(sprintf('K: %d\t rissanen: %f', mtr.K, mtr.rissanen));
   end
   mixture(mtr.K-max(1,finalK)+1) = mtr;
end

if finalK>0
   optmixture = mixture(1);
else
   minriss = mixture(length(mixture)).rissanen; optl=length(mixture);
   for l=length(mixture)-1:-1:1
      if mixture(l).rissanen < minriss
         minriss = mixture(l).rissanen;
         optl = l;
      end
   end
   optmixture = mixture(optl);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PEH EDIT FROM HERE ON

mu = [optmixture.cluster.mu];  % mean of each cluster
in = [find(mu == (min(mu))) find(mu == (max(mu)))]; % index for Gaussian
    % with smaller mean first, then Gaussian with larger mean.
    if length(in) > 2
        in = in(1:2);
    end
mu = mu(in);
[prop_n(1) prop_n(2)] = optmixture.cluster.pb;
prop_n = prop_n(in);

sig = sqrt([optmixture.cluster.R]);    % stdev of each cluster
sig = sig(in);
ll = [optmixture.loglikelihood];

end % GaussianMixture

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
%  Verbatim from original Bowman CLUSTER .m files distributed with
%  CharAnalysis v1.1. Moved here to reduce the distribution from 9 files
%  to 1. No logic changes of any kind.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mixture1 = ClusterNormalize(mixture)

   s = sum([mixture.cluster(:).pb]);
   for k=1:mixture.K
      mixture.cluster(k).pb = mixture.cluster(k).pb/s;
      mixture.cluster(k).invR = inv(mixture.cluster(k).R);
      mixture.cluster(k).const = -(mixture.M*log(2*pi) +log(det(mixture.cluster(k).R)))/2;
   end
   mixture1=mixture;

end

% ?????????????????????????????????????????????????????????????????????????
function mixture1 = EMIterate(mixture, pixels)
% mixture1 = EMIterate(mixture, pixels)
%     perform the EM algorithm with a preassigned fixed order K

[N M] = size(pixels);
Lc = 1+M+0.5*M*(M+1);
epsilon = 0.01*Lc*log(N*M);
[mixture, llnew] = EStep(mixture, pixels);
while true
   llold = llnew;
   mixture = MStep(mixture, pixels);
   [mixture, llnew] = EStep(mixture, pixels);
   if (llnew-llold)<=epsilon
      break;
   end
end
mixture = rmfield(mixture, 'pnk');
mixture.rissanen = -llnew+0.5*(mixture.K*Lc-1)*log(N*M);
mixture.loglikelihood = llnew;
mixture1 = mixture;

end

% ?????????????????????????????????????????????????????????????????????????
function [mixture1 likelihood] = EStep(mixture, pixels)
% EStep     perform the E-step of the EM algorithm
%        1) calculate pnk = Prob(Xn=k|Yn=yn, theta)
%        2) calculate likelihood = log ( prob(Y=y|theta) )

   [N M] = size(pixels);
   K=mixture.K;
   pnk=zeros(N,K);
   for k=1:K
      Y1=pixels-ones(N,1)*mixture.cluster(k).mu';
      Y2=-0.5*Y1*mixture.cluster(k).invR;
      pnk(:,k) = dot(Y1,Y2,2)+mixture.cluster(k).const;
   end
   llmax=max(pnk,[],2);
   pnk =exp( pnk-llmax*ones(1,K) );
   pnk = pnk.*(ones(N,1)*[mixture.cluster(:).pb]);
   ss = sum(pnk,2);
   likelihood = sum(log(ss)+llmax);
   pnk = pnk./(ss*ones(1,K));
   mixture1 = mixture;
   mixture1.pnk = pnk;

end

% ?????????????????????????????????????????????????????????????????????????
function ll = GMClassLikelihood(mixture, Y)
% ll = GMClassLikelihood(class, Y)
% GMClassLikelihood  calculate the log-likelihood of data vectors assuming
% they are generated by the given Gaussian Mixture

   [N, M] = size(Y);
   K = mixture.K;
   pnk=zeros(N,K);
   for k=1:K
      Y1=Y-ones(N,1)*mixture.cluster(k).mu';
      Y2=-0.5*Y1*mixture.cluster(k).invR;
      pnk(:,k) = dot(Y1,Y2,2)+mixture.cluster(k).const;
   end
   llmax=max(pnk,[],2);
   pnk =exp( pnk-llmax*ones(1,K) );
   pnk = pnk.*(ones(N,1)*[mixture.cluster(:).pb]);
   ss = sum(pnk,2);
   ll = log(ss)+llmax;

end

% ?????????????????????????????????????????????????????????????????????????
function mixture = initMixture(pixels, K)

[N M] = size(pixels);
mixture.K=K;
mixture.M=M;
R = (N-1)*cov(pixels)/N;
mixture.Rmin = mean(diag(R))/1e5;
cluster(1).N=0;
cluster(1).pb = 1/K;
cluster(1).mu = pixels(1,:)';
cluster(1).R=R+mixture.Rmin*eye(mixture.M);
if K>1
   period = (N-1)/(K-1);
   for k=2:K
      cluster(k).N=0;
      cluster(k).pb = 1/K;
      cluster(k).mu = pixels(floor((k-1)*period+1),:)';
      cluster(k).R=R+mixture.Rmin*eye(mixture.M);
   end
end
mixture.cluster = cluster;
mixture = ClusterNormalize(mixture);

end

% ?????????????????????????????????????????????????????????????????????????
function mixture1 = MDLReduceOrder(mixture, verbose)
% MDLReduceOrder     reduce the order of the mixture by 1 by combining the
%     two closest distance clusters

   K=mixture.K;
   for k1=1:K
      for k2=k1+1:K
         dist = distance(mixture.cluster(k1), mixture.cluster(k2));
         if (k1==1 && k2==2) || (dist < mindist)
            mink1=k1; mink2=k2;
            mindist = dist;
         end
      end
   end
   if verbose
      disp(['combining cluster: ' int2str(mink1) ' and ' int2str(mink2)]);
   end
   cluster = mixture.cluster;
   cluster(mink1) = addCluster(cluster(mink1), cluster(mink2));
   cluster(mink2:(K-1)) = cluster((mink2+1):K);
   cluster = cluster(1:(K-1));
   mixture1 = mixture;
   mixture1.cluster = cluster;
   mixture1.K = K-1;
   mixture1 = ClusterNormalize(mixture1);

   function d = distance(cluster1, cluster2)
      cluster3 = addCluster(cluster1, cluster2);
      d = cluster1.N*cluster1.const + cluster2.N*cluster2.const-cluster3.N*cluster3.const;
   end

   function cluster3 = addCluster(cluster1, cluster2)
      wt1 = cluster1.N/(cluster1.N+cluster2.N);
      wt2 = 1-wt1;
      M = length(cluster1.mu);
      cluster3 = cluster1;
      cluster3.mu = wt1*cluster1.mu+wt2*cluster2.mu;
      cluster3.R = wt1*(cluster1.R+(cluster3.mu-cluster1.mu)*(cluster3.mu-cluster1.mu)') ...
                     + wt2*(cluster2.R+(cluster3.mu-cluster2.mu)*(cluster3.mu-cluster2.mu)');
      cluster3.invR = inv(cluster3.R);
      cluster3.pb = cluster1.pb+cluster2.pb;
      cluster3.N = cluster1.N+cluster2.N;
      cluster3.const = -(M*log(2*pi) +log(det(cluster3.R)))/2;
   end

end

% ?????????????????????????????????????????????????????????????????????????
function mixture1 = MStep(mixture, pixels)
% MStep     perform the M-step of the EM algorithm
%        from the pnk calculated in the E-step, update parameters of each cluster

   for k=1:mixture.K
      mixture.cluster(k).N = sum(mixture.pnk(:,k));
      mixture.cluster(k).pb = mixture.cluster(k).N;
      mixture.cluster(k).mu = (pixels' * mixture.pnk(:,k))/mixture.cluster(k).N;
      for r=1:mixture.M
         for s=r:mixture.M
            mixture.cluster(k).R(r,s) = ((pixels(:,r)-mixture.cluster(k).mu(r))' ...
                           * ((pixels(:,s)-mixture.cluster(k).mu(s)).*mixture.pnk(:,k))) ...
                           /mixture.cluster(k).N;
            if r~=s
               mixture.cluster(k).R(s,r) = mixture.cluster(k).R(r,s);
            end
         end
      end
      mixture.cluster(k).R = mixture.cluster(k).R+mixture.Rmin*eye(mixture.M);
   end
   mixture1=ClusterNormalize(mixture);
end