function dh = hellingerPoisson(alpha, beta)
%hellingerPoisson Squared Hellinger distance between two Poisson distributions
%   dh = hellingerPoisson(alpha, beta)
%   
%   Squared Hellinger distance (dh) between two Poisson distributions with
%   rates alpha and beta.

assert(ndims(alpha) == ndims(beta) && all(size(alpha) == size(beta)), 'alpha and beta must be the same size')

dh = 1 - exp(-0.5 .* (sqrt(alpha) - sqrt(beta)).^2);

end

