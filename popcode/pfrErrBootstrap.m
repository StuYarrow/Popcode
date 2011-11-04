function [pfr, pfrSEM] = pfrErrBootstrap(f, r, ns)

% Get means
fMu = cellfun(@(a) mean(a.samples), f, 'UniformOutput', false);
rMu = cellfun(@(a) mean(a.samples), r, 'UniformOutput', false);
mMu = cellfun(@minus, fMu, rMu, 'UniformOutput', false);
pfr = cellfun(@(a) 2 * a(2) / (a(1) + a(3)), mMu);

% Get marginal samples
mSamps = cellfun(@(a,b) a.samples - b.samples, f, r, 'UniformOutput', false);

clear f r

% Compute mean, SEM
mMu = cellfun(@mean, mSamps, 'UniformOutput', false);
mStd = cellfun(@(a) sqrt(var(a) / size(a,1)), mSamps, 'UniformOutput', false);

% Preallocate
samps = zeros(size(mMu,1), size(mMu,2), ns);

for i = 1 : ns
    % Sample
    m = cellfun(@(a,b) a + b .* randn(size(a)), mMu, mStd, 'UniformOutput', false);
    
    % Compute and store PFR sample
    samps(:,:,i) = cellfun(@(a) 2 * a(2) / (a(1) + a(3)), m);
end

% Compute StdErr
pfrSEM = std(samps,0,3);
