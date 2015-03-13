function mLH = getMarginalLikelihoodBGe(A,D,V)
	ns = size(A,1);
	nu = 1;
	alpha = ns + 2;
	mu0 = zeros(1,ns);
	T0 = (nu*(alpha-ns-1)/(nu+1)) * eye(ns);
	
	mLH = zeros(ns,1);
	for i=1:ns
		pidx = find(A(:,i)==1)';
		Dv = D(V(:,i)==1,:);
		Di = Dv(:,i);
		Dpi = Dv(:,pidx);
		mLH(i) = logLH([Di,Dpi],nu,alpha,mu0([i pidx]),T0([i pidx],[i pidx])) - ...
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