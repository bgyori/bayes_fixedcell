function y = smoothkernel(X,Y,x)
	b = 15;
	y = zeros(size(x));
	for i=1:length(x)
		Kx = exp(-(x(i)-X).^2/(2*b^2));
		y(i) = sum(Kx.*Y)/sum(Kx);
	end
end