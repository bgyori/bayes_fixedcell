function x = xfromhist(h,s)
	hn = h/sum(h);
	hs = hn/min(hn(hn>0));
	counts = round(hs);
	x = [];
	for i=1:length(h)
		x = [x repmat(s(i),1,counts(i))];
	end
end