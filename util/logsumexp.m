function s = logsumexp(X)
	if isvector(X)
		X = X(:)';
	end
	m = max(X,[],2);
	Y = X - repmat(m,1,size(X,2));
	s = m + log(sum(exp(Y),2));
end