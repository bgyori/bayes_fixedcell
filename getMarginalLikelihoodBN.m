function mLH = getMarginalLikelihoodBN(A,S,D,V)
% Marginal likelihood for binary Bayesian network
	nLevel = 2;
	ns = size(A,1);
	
	if size(D,2)~=ns || size(V,2)~=ns
		error('D or V has wrong size');
	end
	
	mLH = zeros(1,ns);
	for i=1:ns
		L = 0;
		Dv = D(V(:,i)==1,:);
		nodeData = Dv(:,i);
		
		% Get parents of node
		parentIdx = find(A(:,i));
		% Number of parents
		nParents = length(parentIdx);
		if nParents == 0
			dataSamp = size(Dv,1);
			priorSamp = S/nLevel;
			L = L + gammaln(S) - gammaln(S+dataSamp);
			for k=0:nLevel-1
				dataSamp = sum(nodeData==k);
				L = L + gammaln(priorSamp+dataSamp) - gammaln(priorSamp);
			end
		else
			% Number of configurations
			nConfig = 2^nParents;
			% Prior sample size
			priorSamp = S/(nLevel*nConfig);
			parentData = Dv(:,parentIdx);
			configs = binmat2dec(parentData);
			for j=0:nConfig-1
				configDataRows = (configs==j);
				if any(configDataRows)
					dataSamp = sum(configDataRows);
					L = L + gammaln(priorSamp*nLevel)-gammaln(priorSamp*nLevel+dataSamp);
					for k=0:nLevel-1
						dataSamp = sum(nodeData(configDataRows)==k);
						L = L + gammaln(priorSamp+dataSamp)-gammaln(priorSamp);
					end
				end
			end
		end
		mLH(i) = L;
	end
end