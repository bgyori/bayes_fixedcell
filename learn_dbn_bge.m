function mLH = learn_dbn_bge(A, Dminus, Dplus, V)
	ns = size(A, 1);
	mLH = zeros(1, ns);
	
	nu = 1;
	alpha = ns + 2;
	mu0 = zeros(1,ns);
	sigma = (nu*(alpha-ns-1)/(nu+1));
	T0 = sigma * eye(ns);
	
	for i=1:ns
		Dvplus = Dplus(:,V(i,:)==1);
		Dvminus = Dminus(:,V(i,:)==1);
		parentIdx = find(A(:,i)==1)';
		nodeData = Dvplus(i,:);
		parentData = Dvminus(parentIdx,:);
	
		T0i = T0([i,parentIdx],[i,parentIdx]);
		mu0i = mu0([i parentIdx]);
		mLH(i) = logLH([nodeData; parentData], nu, alpha, mu0i, T0i) - ...
				 logLH(parentData, nu, alpha, mu0(parentIdx), ...
					   T0(parentIdx, parentIdx));
		priorG = -log(nchoosek(ns-1,length(parentIdx)));
		mLH(i) = mLH(i) + priorG;
	end
end

function L = logLH(D, nu, alpha, mu0, T0)
	if isempty(D)
		L = 0;
		return;
	end
	[N,m] = size(D);
	M = mean(D');
	S = (m-1)*cov(D');
	T = T0 + S + nu*m/(nu+m)*(mu0 - M)*(mu0 - M)';
	% Terms of the log-likelihood
	L1 = -0.5 * N * m * log(2*pi);
	L2 = 0.5 * N * (log(nu) - log(nu+m));
	L3 = logC(N, alpha);
	L4 = -logC(N, alpha+m);
	L5 = 0.5 * alpha * log(det(T0));
	L6 = -0.5 * (alpha + m) * log(det(T));
	% Add up all the terms
	L = L1 + L2 + L3 + L4 + L5 + L6;
end

function C = logC(N, alpha)
	% Terms of log(C)
	C1 = -0.5 * alpha * N * log(2);
	C2 = -0.25 * N * (N-1) * log(pi);
	C3 = -sum(gammaln(0.5 * (alpha + 1 - (1:N))));
	C = C1 + C2 + C3;
end