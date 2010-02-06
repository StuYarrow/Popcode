function r = meanR(obj, stim)
% CIRCGAUSSNEURONS/MEANR calculates mean responses to a set of stimuli
% r = meanR(obj, stimulusEnsemble)
%
% r = maxRate * exp(-((stimulus - preferredStimulus)^2 / (2 * width^2))) + backgroundRate

if ~isa(stim, 'StimulusEnsemble') 
	error([inputname(2) ' is not a valid StimulusEnsemble object'])
end

if stim.dimensionality ~= 1
	error('CircGaussNeurons only supports 1-D stimuli at present')
end

stims = repmat(stim.ensemble, obj.Neurons.popSize, 1);

maxRate = repmat(obj.Neurons.maxRate, 1, stim.n);
backgroundRate = repmat(obj.Neurons.backgroundRate, 1, stim.n);
centre = repmat(obj.Neurons.preferredStimulus, 1, stim.n);
width = repmat(obj.width, 1, stim.n);

r = maxRate .* exp(-(1.0 ./ deg2rad(width).^2) .* (1 - cosd(double(stims - centre)))) + backgroundRate;