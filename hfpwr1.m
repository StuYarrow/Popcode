function hfp = hfpwr1(sig)
%
%
%

S = abs(fft(sig, [], 2));
S = S ./ repmat(max(S, [], 2), [1 size(S, 2)]);
S = S(:,floor(end/4):floor(3*end/4));
hfp = mean(S, 2);