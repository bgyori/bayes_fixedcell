function [mLH,E] = learn_dbn_gauss(A,Dminus,Dplus,V,int_terms)
	ns = size(A,1);
 	E = zeros(ns);
	
	mLH = zeros(1, ns);
	for i=1:ns
		Dvminus = Dminus(:,V(i,:)==1);
		Dvplus = Dplus(:,V(i,:)==1);
		X = Dvplus(i,:)';
		X = stdize(X);
		
		if isempty(Dvminus)
			mLH(i) = nan;
			continue;
		end
		
		n = size(Dvplus,2);
		
		parents = find(A(:,i));
		np = length(parents);
		
		
		if int_terms
			npc = 2^np;
			c = (1+n)^(-(npc-1)/2);
			B = zeros(n,npc-1);
			B1 = Dvminus(parents,:)';
			for k=1:npc-1
				mask = dec2binvec(k,np);
				B(:,k) = prod(B1(:,mask),2);
			end
		else
			c = (1+n)^(-np/2);
			B = Dvminus(parents,:)';
		end
		
		B = stdize(B);
		
		BB = B'*B;
		if cond(BB) > 1e4
			BB = BB + 0.1*eye(size(BB));
		end
		Bpinv = (B*inv(BB))*B';
		
		
		mLH(i) = log(c) + (-n/2)*log((X'*X - n/(n+1)*X'*Bpinv*X));
		
		
		for j=1:ns
 			cc = corrcoef(Dvminus(j,:),X);
 			E(j,i) = cc(1,2);
 		end
	end
end

function Y = stdize(Y)
	Y = Y-repmat(mean(Y),size(Y,1),1);
	Y = Y./max(repmat(std(Y),size(Y,1),1),1e-10);
end