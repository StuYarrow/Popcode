function [pfr, pfrSEM] = pfrErr(f, r)

fPeak = cellfun(@(a) a.runMean(2), f);
fPeakSEM = cellfun(@(a) a.runSEM(2), f);

rPeak = cellfun(@(a) a.runMean(2), r);
rPeakSEM = cellfun(@(a) a.runSEM(2), r);

mPeak = fPeak - rPeak;
mPeakSEM = sqrt(fPeakSEM.^2 + rPeakSEM.^2);

fFlank = cellfun(@(a) a.runMean([1 3]), f, 'UniformOutput', false);
fFlankSEM = cellfun(@(a) a.runSEM([1 3]), f, 'UniformOutput', false);
fFlank = cellfun(@mean, fFlank, repmat({2}, size(fFlank)));
fFlankSEM = cellfun(@(a) mean(a,2) / sqrt(2), fFlankSEM);

rFlank = cellfun(@(a) a.runMean([1 3]), r, 'UniformOutput', false);
rFlankSEM = cellfun(@(a) a.runSEM([1 3]), r, 'UniformOutput', false);
rFlank = cellfun(@mean, rFlank, repmat({2}, size(rFlank)));
rFlankSEM = cellfun(@(a) mean(a,2) / sqrt(2), rFlankSEM);

mFlank = fFlank - rFlank;
mFlankSEM = sqrt(fFlankSEM.^2 + rFlankSEM.^2);

pfr = mPeak ./ mFlank;
pfrSEM = pfr .* sqrt((mPeakSEM ./ mPeak).^2 + (mFlankSEM ./ mFlank).^2);
