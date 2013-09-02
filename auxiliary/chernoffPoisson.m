function dc = chernoffPoisson(alpha, beta)
%chernoffPoisson Chernoff distance between two Poisson distributions
%   dh = chernoffPoisson(alpha, beta)
%   
%   Chernoff distance (dh) between two Poisson distributions with
%   rates alpha and beta.
%   
%   After: D. H. Johnson and S. Sinanovic, "Symmetrizing the Kullback-Leibler distance",  2001.

assert(ndims(alpha) == ndims(beta) && all(size(alpha) == size(beta)), 'alpha and beta must be the same size')

r = beta ./ alpha;
dc = alpha .* ((r - 1) .* (log((r - 1) ./ log(r)) - 1) + log(r)) ./ log(r);

end

