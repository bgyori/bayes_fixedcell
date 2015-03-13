function mLH = getMarginalLikelihoodDBGe(A,D,V)
	[ns,nt,nc] = size(D);
	nu = 1;
	alpha = ns + 2;
	mu0 = zeros(1,ns);
	sigma = (nu*(alpha-ns-1)/(nu+1));
	T0 = sigma * eye(ns);

	
	mLH = zeros(ns,1);
	for i=1:ns
		pidx = find(A(:,i)==1)';
		Dv = D(:,:,V(i,:)==1);
		Di = reshape(Dv(i,2:end,:),[1,(nt-1)*size(Dv,3)])';
		Dpi = reshape(Dv(:,1:end-1,:),[ns,(nt-1)*size(Dv,3)]);
		Dpi = Dpi(pidx,:)';
		%Modify for the case when node is its own parent
		if any(pidx==i)
			T0i =  sigma * eye(length(pidx)+1);
			mu0i = zeros(1,length(pidx)+1);
		else
			T0i = T0([i pidx],[i pidx]);
			mu0i = mu0([i pidx]);
		end
		mLH(i) = logLH([Di,Dpi],nu,alpha,mu0i,T0i) - ...
				 logLH(Dpi,nu,alpha,mu0(pidx),T0(pidx,pidx));
	end
end

function L = logLH(D,nu,alpha,mu0,T0)
	if isempty(D)
		L = 0;
		return;
	end
	[m,N] = size(D);
	M = mean(D);
	S = (D-repmat(M,m,1))'*(D-repmat(M,m,1));
	T = T0 + S + nu*m/(nu+m)*(mu0 - M)*(mu0 - M)';
	L = -0.5*N*m*log(2*pi) + 0.5*N*(nu/(nu+m));
	L = L + logC(N,alpha) - logC(N,alpha+m);
	L = L + 0.5*alpha*log(det(T0)) - 0.5*(alpha+m)*log(det(T));
end

function C = logC(N,a)
	C = 0.5*a*N*log(2) + 0.25*N*(N-1)*log(pi);
	C = C + sum(gammaln(0.5*(a+1-(1:N))));
end