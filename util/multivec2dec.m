function dec = multivec2dec(x, levels)
% MULTIVEC2DEC: Converts a multinomial vector integer into a decimal integer
% X: Multinomial vector eg. [1 0 3], least->most significant bit
% LEVELS: Number of levels of each position in the multinomial vector e.g.
% [2 1 4]
% DEC: Decimal integer e.g. 7
	cl = [1, cumprod(levels(1:end-1))];
	dec = sum(x(:) .* cl(:));
end