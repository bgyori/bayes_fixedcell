function mLH = learn_dbn_bde(A,Dminus,Dplus,V,S,levels)
% Marginal likelihood for binary dynamic Bayesian network
	ns = size(A,1);
	mLH = zeros(1,ns);
	for i=1:ns
        Dvplus = Dplus(:,V(i,:)==1);
		Dvminus = Dminus(:,V(i,:)==1);
		parentIdx = find(A(:,i)==1);
		nodeData = Dvplus(i,:);
		parentData = Dvminus(parentIdx,:);
		nlevels = levels(i);
		plevels = levels(parentIdx);
		mLH(i) = bde_counts_multi(parentData, nodeData, plevels, nlevels, S);
	end
end