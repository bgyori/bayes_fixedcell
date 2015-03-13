function [mLH,E] = getMarginalLikelihood(A,D,V)
	intTerms = true;

	[ns,nt,nc] = size(D);
	
	if nargin < 3
		V = ones(ns,nc);
	end
	
 	E = zeros(ns);
	
	
	mLH = zeros(1,ns);
	for i=1:ns
		Dv = D(:,:,V(i,:)==1);
		X = reshape(Dv(i,2:end,:),[1,(nt-1)*size(Dv,3)])';
		X = stdize(X);
		Dminus = reshape(Dv(:,1:end-1,:),[ns,(nt-1)*size(Dv,3)]);
		
		if isempty(Dminus)
			mLH(i) = nan;
			continue;
		end
		
		n = size(Dv,3)*(nt-1);
		
		parents = find(A(:,i));
		np = length(parents);
		
		
		if intTerms
			npc = 2^np;
			c = (1+n)^(-(npc-1)/2);
			B = zeros(n,npc-1);
			B1 = Dminus(parents,:)';
			for k=1:npc-1
				mask = dec2binvec(k,np);
				B(:,k) = prod(B1(:,mask),2);
			end
		else
			c = (1+n)^(-np/2);
			B = Dminus(parents,:)';
		end
		
		B = stdize(B);
		
		BB = B'*B;
		if cond(BB) > 1e4
			BB = BB + 0.1*eye(size(BB));
		end
		Bpinv = (B*inv(BB))*B';
		
		
		mLH(i) = log(c) + (-n/2)*log((X'*X - n/(n+1)*X'*Bpinv*X));
		
		
		for j=1:ns
 			cc = corrcoef(Dminus(j,:),X);
 			E(j,i) = cc(1,2);
 		end
	end
end

function Y = stdize(Y)
	Y = Y-repmat(mean(Y),size(Y,1),1);
	Y = Y./max(repmat(std(Y),size(Y,1),1),1e-10);
end