function gn = CircGaussNeurons(varargin)
% CIRCGAUSSNEURONS/CIRCGAUSSNEURONS Constructor for CircGaussNeurons object - population of neurons with circular Gaussian tuning curves
% n = CircGaussNeurons(popSize, preferredStimulus, width, maxRate, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
%
% popSize - population size (number of neurons)
% dimensionality - stimulus dimensionality (only 1-D stimuli currently supported)
% preferredStimulus - preferred stimulus value
% width - tuning curve width (this would be the variance if the curve was a probability distribution)
% maxRate - maximum firing rate (Hz)
% backgroundRate - background (spontaneous) firing rate (Hz)
% integrationTime - spike counting time per trial
% variability - a Variability object
% variabilityOpts - vector of options
%
% preferredStimulus, width, maxFiringRate and backgroundFiringRate can be scalars or vectors of length popSize.
% CircGaussNeurons accept only 1-D stimuli at present.

switch nargin
case 0
	% Create default object
	n = Neurons();
	gn.width = [];

case 1
	% If passed a valid CircGaussNeurons object, return it
	if isa(in1, 'CircGaussNeurons')
		gn = in1;
	else
		error([inputname(1) ' is not a valid CircGaussNeurons object'])
	end

case 7
	% Default constructor
	n = Neurons(1, varargin{1}, varargin{3}, varargin{4}, varargin{5}, varargin{6}, varargin{7});
	
	if length(varargin{2}) == 1 & isnumeric(varargin{2})
		gn.width = double(varargin{2}(ones(n.popSize, 1)));
	elseif length(varargin{2}) == n.popSize & isvector(varargin{2}) & isnumeric(varargin{2})
		gn.width = reshape(double(), n.popSize, 1);
	else
		error([inputname(3) ' is not a valid maximum firing rate value or vector for population size ' n.popSize])
	end

otherwise
	error('Wrong number of arguments')
end

gn = class(gn, 'CircGaussNeurons', n);