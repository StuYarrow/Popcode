function plt = poplot(stim, measure)
%
%
%

stims = double(stim.ensemble);
[dummy, i] = sort(stims);
lines = plot(stims(i), measure(i));

plt = get(lines, 'Parent');
set(plt, 'XLim', [min(stims) max(stims)])
