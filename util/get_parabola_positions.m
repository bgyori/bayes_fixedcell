function pos = get_parabola_positions(p, x, y)
	pos = zeros(length(x), 1);
	for i=1:length(x)
		b = [x(i), y(i)];
		fobj = @(x) get_parabola_distance(p, x, b);
		posx = fminsearch(fobj, 0.5);
		pos(i) = distance_along_parabola(p, posx);
	end
end

function pos = distance_along_parabola(p, x)
	% p is a vector containing the parameters of the parabola
	% x is a point on the x-axis
	x1 = (-p(2)-sqrt(p(2)*p(2) - 4*p(1)*p(3)))/(2*p(1));
	x2 = (-p(2)+sqrt(p(2)*p(2) - 4*p(1)*p(3)))/(2*p(1));
	xmin = min([x1, x2]);
	xmax = max([x1, x2]);
	total_arc = get_parabola_arc_length(p, xmin, xmax);
	arc = get_parabola_arc_length(p, xmin, x);
	pos = arc/total_arc;
end

function d = get_parabola_arc_length(p, x1, x2)
	fun = @(x) sqrt(1 + (p * [x.^2; x; ones(1, length(x))]).^2);
	d = integral(fun, x1, x2);
end
	
function d = get_parabola_distance(p, x, b)
	% p is a vector containing the parameters of the parabola
	% X is a point on the x-axis
	% b is a point in space that we want to map
	
	% Point on the parabola at position x
	a = [x, dot(p, [x^2, x, 1])];
	% Length of vector from a to b
	d = norm(a-b);
end