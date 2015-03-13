function [h,mu,sigma] = emgauss(x,m)
	x = x(:);
	n = length(x);
	h = 0.5*ones(n,m);
	mu = linspace(0.2,0.8,m);
	sigma = 0.1*ones(m,1);
	
	t = 1;
	while t < 10

		% Maximization
		for i=1:m
			w(i) = sum(h(:,i))/n;
		end
		
		for j=1:n
			for i=1:m
				h(j,i) = w(i)*normpdf(x(j),mu(i),sigma(i));
			end
			h(j,:) = h(j,:)/sum(h(j,:));
		end
		t = t + 1;
		
		% Expectation
		for i=1:m
			
			mu(i) = sum(h(:,i).*x)/sum(h(:,i));
			sigma(i) = sqrt(sum(h(:,i).*(x-mu(i)).*(x-mu(i)))/sum(h(:,i)));
		end

	end
end