function signalNorm = normalize_foxo(signal)
	% Make row vectors
	medians = signal.qm(:)';
	iqrs = signal.iqr(:)';
	% Sort medians in ascending order and sort the iqrs accordingly
	[shiftmedians, indsorted] = sort(medians);
	shiftiqrs = iqrs(indsorted);
	
	% Take the highest 8 values of the medians and
	% corresponding iqrs
	igfhighmedians = shiftmedians(end-7:end);
	igfhighiqrs = shiftiqrs(end-7:end);
	% Take the lowest 8 values of the medians and
	% the corresponding iqrs
	shiftmedians = shiftmedians(1:8);
	shiftiqrs = shiftiqrs(1:8);
	% Take the median of the selected medians and
	% the median of the selected iqrs. This becomes
	% the shift value. 
	shiftpar = median(shiftmedians);
	shiftiqr = median(shiftiqrs);
	% Deduct the median shift and the iqr shift from the row vectors
	medians = medians - shiftpar;
	iqrs = iqrs - shiftiqr;
	% Deduct the median shift and iqr from the highest medians and 
	% corresponding iqrs
	igfhighmedians = igfhighmedians - shiftpar;
	igfhighiqrs = igfhighiqrs - shiftiqr;

	% Apply the same to the original matrix
	signalNorm.qm = signal.qm - shiftpar;
	signalNorm.iqr = signal.iqr - shiftiqr;
	
	% Get the angle with the median of the highest medians along X axis
	% and median of the corresponding iqrs along the Y axis
	alpha = atan(median(igfhighiqrs)/median(igfhighmedians));
	% Define the rotation matrix for rotating clockwise by alpha
	Rmat = [cos(alpha)  sin(alpha); ...
			-sin(alpha) cos(alpha)];
	% Rotate both medians and iqrs
	tmp = Rmat*[medians; iqrs];
	% Reconstruct the now rotated row matrices
	medians = tmp(1,:);
	iqrs = tmp(2,:);

	% Do the same rotation for the original matrix
	s = size(signalNorm.qm);
	tmp2 = Rmat*[signalNorm.qm(:)';signalNorm.iqr(:)'];
	signalNorm.qm = reshape(tmp2(1,:)',s);
	signalNorm.iqr = reshape(tmp2(2,:)',s);

	% Take the highest medians again, after rotation
	igfhighmedians = sort(medians);
	igfhighmedians = igfhighmedians(end-7:end);

	% Divide the medians after rotation by the median of the
	% highest medians
	signalNorm.qm = signalNorm.qm/median(igfhighmedians);
end