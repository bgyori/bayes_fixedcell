function y = auc(x,t)
	nc = size(x,2);
	dt = t(2:end)-t(1:end-1);
	dt = repmat(dt(:),1,nc);
	y = sum(0.5*dt.*(x(1:end-1,:)+x(2:end,:)));
end