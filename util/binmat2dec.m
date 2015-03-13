function dec = binmat2dec(B)
	[r,c] = size(B);
	dec = sum(B.*repmat(2.^(0:c-1),r,1),2);
end