function pos = get_parabola_positions(p, x, y)
	pos = zeros(length(x), 1);
	for i=1:length(x)
		b = [x(i), y(i)];
		fobj = @(x) get_parabola_distance(x, b, p);
		pos(i) = fminsearch(fobj, 0.5);
	end
end

function d = get_parabola_distance(x, b, p)
	% X is a point on the x-axis
	% b is a point in space that we want to map
	% p is a vector containing the parameters of the parabola
	
	% Point on the parabola at position x
	a = [x, dot(p, [x^2, x, 1])];
	% Length of vector from a to b
	d = norm(a-b);
end