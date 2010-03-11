% Test cases script

zap
orient portrait
tol = 1e-3;

%dctconfig('hostname', '127.0.0.1');
%warning('off', 'MATLAB:nearlySingularMatrix')

stim = StimulusEnsemble('circular', 360, 180);
ds = datestr(today, 'yyyymmdd');

% Test Case 1 (SSI)
% "Effect of varying beta" EC Dissertation, p 30
tc1a = CircGaussNeurons([0]', 30.0, 1.0, 0.0, 1.0, 'Gaussian-independent', [0.026 0.5 0.024]);
tc1b = CircGaussNeurons([0]', 30.0, 1.0, 0.0, 1.0, 'Gaussian-independent', [0.21 0.5 0.024]);
tc1c = CircGaussNeurons([0]', 30.0, 1.0, 0.0, 1.0, 'Gaussian-independent', [0.42 0.5 0.024]);
ssi1a = ssi(tc1a, [], 'randMC', stim, tol, 4000);
ssi1b = ssi(tc1b, [], 'randMC', stim, tol, 4000);
ssi1c = ssi(tc1c, [], 'randMC', stim, tol, 4000);
clf
plot(double(stim.ensemble), [ssi1a' ssi1b' ssi1c'] ./ max(max([ssi1a' ssi1b' ssi1c'])))
axis([-180 180 0 1.4])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 1')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
orient portrait
print('-depsc', [ds '_tc1.eps'])


% Test Case 2 (SSI)
% PS - "Fisher vs Shannon...", AREADNE June 2008, fig 1
tc2a = CircGaussNeurons([0]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [0.2 0.5]);
tc2b = CircGaussNeurons([0]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [2.0 0.5]);
fi2a = fisher(tc2a, 'analytic', stim, 0);
ssi2a = ssi(tc2a, [], 'randMC', stim, tol, 4000);
fi2b = fisher(tc2b, 'analytic', stim, 0);
ssi2b = ssi(tc2b, [], 'randMC', stim, tol, 4000);
clf
plot(subplot(1,2,1), double(stim.ensemble), ssi2a ./ max(ssi2a), 'g-', double(stim.ensemble), fi2a ./ max(fi2a), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 2 F=0.2')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,2,2), double(stim.ensemble), ssi2b ./ max(ssi2b), 'g-', double(stim.ensemble), fi2b ./ max(fi2b), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 2 F=2.0')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
orient portrait
print('-depsc', [ds '_tc2.eps'])


% Test Case 3 (Marginal SSI)
% PS - "Fisher vs Shannon...", AREADNE June 2008, fig 2
tc3fi = CircGaussNeurons([0]', 33.0, 50.0, 0.0, 0.015, 'Gaussian-independent', [1.0 0.5]);
tc3a = CircGaussNeurons([-180 -90 0 90]', 33.0, 50.0, 0.0, 0.015, 'Gaussian-independent', [1.0 0.5]);
tc3b = CircGaussNeurons([-180 : 22.5 : 180-22.5]', 33.0, 50.0, 0.0, 0.015, 'Gaussian-independent', [1.0 0.5]);
tc3c = CircGaussNeurons([-180 : 12 : 180-12]', 33.0, 50.0, 0.0, 0.015, 'Gaussian-independent', [1.0 0.5]);
fi3 = fisher(tc3fi, 'analytic', stim, 0);
ssi3a = ssi(tc3a, [3], 'randMC', stim, tol, 4000);
ssi3b = ssi(tc3b, [9], 'randMC', stim, tol, 10000);
ssi3c = ssi(tc3c, [16], 'randMC', stim, tol, 20000);
clf
plot(subplot(1,3,1), double(stim.ensemble), ssi3a ./ max(ssi3a), 'g-', double(stim.ensemble), fi3 ./ max(fi3), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 3 N=4')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,3,2), double(stim.ensemble), ssi3b ./ max(ssi3b), 'g-', double(stim.ensemble), fi3 ./ max(fi3), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 3 N=16')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,3,3), double(stim.ensemble), ssi3c ./ max(ssi3c), 'g-', double(stim.ensemble), fi3 ./ max(fi3), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 3 N=30')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
print('-depsc', [ds '_tc3.eps'])


% Test Case 4 (Marginal SSI/Fisher)
% PS - AREADNE 2008 Poster, Fig 3
tc4fia = CircGaussNeurons([0]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [1.0 0.5]);
tc4fib = CircGaussNeurons([-36]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [1.0 0.5]);
tc4fid = CircGaussNeurons([-25.7]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [1.0 0.5]);
tc4a = CircGaussNeurons([-180 -90 0 90]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [1.0 0.5]);
tc4b = CircGaussNeurons([-180 : 360/5 : 180-360/5]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [1.0 0.5]);
tc4c = CircGaussNeurons([-180 : 360/6 : 180-360/6]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [1.0 0.5]);
tc4d = CircGaussNeurons([-180 : 360/7 : 180-360/7]', 33.0, 50.0, 0.0, 1.0, 'Gaussian-independent', [1.0 0.5]);
fi4a = fisher(tc4fia, 'analytic', stim, 0);
fi4b = fisher(tc4fib, 'analytic', stim, 0);
fi4d = fisher(tc4fid, 'analytic', stim, 0);
ssi4a = ssi(tc4a, [3], 'randMC', stim, tol, 10000);
ssi4b = ssi(tc4b, [3], 'randMC', stim, tol, 10000);
ssi4c = ssi(tc4c, [4], 'randMC', stim, tol, 10000);
ssi4d = ssi(tc4d, [4], 'randMC', stim, tol, 10000);
clf
plot(subplot(1,4,1), double(stim.ensemble), ssi4a ./ max(ssi4a), 'g-', double(stim.ensemble), fi4a ./ max(fi4a), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 4 N=4')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,2), double(stim.ensemble), ssi4b ./ max(ssi4b), 'g-', double(stim.ensemble), fi4b ./ max(fi4b), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 4 N=5')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,3), double(stim.ensemble), ssi4c ./ max(ssi4c), 'g-', double(stim.ensemble), fi4a ./ max(fi4a), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 4 N=6')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,4), double(stim.ensemble), ssi4d ./ max(ssi4d), 'g-', double(stim.ensemble), fi4d ./ max(fi4d), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 4 N=7')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
print('-depsc', [ds '_tc4.eps'])


% Test Case 5 (Marginal SSI/Fisher)
% PS - AREADNE 2008 Poster, Fig 4
tc5fi = CircGaussNeurons([0]', 33.0, 50.0, 0.0, 0.1, 'Gaussian-independent', [1.0 0.5]);
tc5a = CircGaussNeurons([-180 -90 0 90]', 33.0, 50.0, 0.0, 0.1, 'Gaussian-independent', [1.0 0.5]);
tc5b = CircGaussNeurons([-180 : 360/8 : 180-360/8]', 33.0, 50.0, 0.0, 0.1, 'Gaussian-independent', [1.0 0.5]);
tc5c = CircGaussNeurons([-180 : 360/16 : 180-360/16]', 33.0, 50.0, 0.0, 0.1, 'Gaussian-independent', [1.0 0.5]);
tc5d = CircGaussNeurons([-180 : 360/30 : 180-360/30]', 33.0, 50.0, 0.0, 0.1, 'Gaussian-independent', [1.0 0.5]);
fi5 = fisher(tc5fi, 'analytic', stim, 0);
ssi5a = ssi(tc5a, [3], 'randMC', stim, tol, 10000);
ssi5b = ssi(tc5b, [5], 'randMC', stim, tol, 10000);
ssi5c = ssi(tc5c, [9], 'randMC', stim, tol, 40000);
ssi5d = ssi(tc5d, [16], 'randMC', stim, tol, 40000);
clf
plot(subplot(1,4,1), double(stim.ensemble), ssi5a ./ max(ssi5a), 'g-', double(stim.ensemble), fi5 ./ max(fi5), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 5 N=4')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,2), double(stim.ensemble), ssi5b ./ max(ssi5b), 'g-', double(stim.ensemble), fi5 ./ max(fi5), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 5 N=8')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,3), double(stim.ensemble), ssi5c ./ max(ssi5c), 'g-', double(stim.ensemble), fi5 ./ max(fi5), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 5 N=16')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,4), double(stim.ensemble), ssi5d ./ max(ssi5d), 'g-', double(stim.ensemble), fi5 ./ max(fi5), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 5 N=30')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
print('-depsc', [ds '_tc5.eps'])


% Test Case 6 (Marginal SSI/Fisher)
% PS - AREADNE 2008 Poster, Fig 5
tc6fi = CircGaussNeurons([0]', 33.0, 50.0, 0.0, 0.02, 'Gaussian-independent', [1.0 0.5]);
tc6a = CircGaussNeurons([-180 -90 0 90]', 33.0, 50.0, 0.0, 0.02, 'Gaussian-independent', [1.0 0.5]);
tc6b = CircGaussNeurons([-180 : 360/8 : 180-360/8]', 33.0, 50.0, 0.0, 0.02, 'Gaussian-independent', [1.0 0.5]);
tc6c = CircGaussNeurons([-180 : 360/16 : 180-360/16]', 33.0, 50.0, 0.0, 0.02, 'Gaussian-independent', [1.0 0.5]);
tc6d = CircGaussNeurons([-180 : 360/30 : 180-360/30]', 33.0, 50.0, 0.0, 0.02, 'Gaussian-independent', [1.0 0.5]);
fi6 = fisher(tc6fi, 'analytic', stim, 0);
ssi6a = ssi(tc6a, [3], 'randMC', stim, tol, 10000);
ssi6b = ssi(tc6b, [5], 'randMC', stim, tol, 10000);
ssi6c = ssi(tc6c, [9], 'randMC', stim, tol, 40000);
ssi6d = ssi(tc6d, [16], 'randMC', stim, tol, 40000);
clf
plot(subplot(1,4,1), double(stim.ensemble), ssi6a ./ max(ssi6a), 'g-', double(stim.ensemble), fi6 ./ max(fi6), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 6 N=4')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,2), double(stim.ensemble), ssi6b ./ max(ssi6b), 'g-', double(stim.ensemble), fi6 ./ max(fi6), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 6 N=8')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,3), double(stim.ensemble), ssi6c ./ max(ssi6c), 'g-', double(stim.ensemble), fi6 ./ max(fi6), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 6 N=16')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
plot(subplot(1,4,4), double(stim.ensemble), ssi6d ./ max(ssi6d), 'g-', double(stim.ensemble), fi6 ./ max(fi6), 'b--')
axis([-180 180 0 1.2])
set(gca, 'XTick', [-180 -90 0 90 180])
title('Test Case 6 N=30')
xlabel('Stimulus (deg)')
ylabel('SSI (bits)')
print('-depsc', [ds '_tc6.eps'])
