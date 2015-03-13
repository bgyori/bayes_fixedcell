function  [thresh,em] = otsu(x)
  num_bins = 20;
  bins = linspace(min(x),max(x),num_bins);
  counts = hist(x,num_bins);
  
  % Variables names are chosen to be similar to the formulas in
  % the Otsu paper.
  p = counts / sum(counts);
  omega = cumsum(p);
  mu = cumsum(p .* bins);
  mu_t = mu(end);
  
  sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));

  % Find the location of the maximum value of sigma_b_squared.
  % The maximum may extend over several bins, so average together the
  % locations.  If maxval is NaN, meaning that sigma_b_squared is all NaN,
  % then return 0.
  maxval = max(sigma_b_squared);
  isfinite_maxval = isfinite(maxval);
  if isfinite_maxval
    thresh = mean(bins(sigma_b_squared == maxval));
    % Normalize the threshold to the range [0, 1].
    %thresh = (idx - 1) / (num_bins - 1);
  else
    thresh = 0.0;
  end

  % compute the effectiveness metric
  em = maxval/(sum(p.*(bins.^2)) - mu_t^2);
end