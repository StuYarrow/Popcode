function varargout = get(obj, prop_name)
% GET Get Neurons properties from the specified object and return the value

switch prop_name
case 'Neurons'
	varargout = {obj.Neurons};
case 'popSize'
	varargout = {obj.Neurons.popSize};
case 'dimensionality'
    varargout = {obj.Neurons.dimensionality};
case 'preferredStimulus'
	varargout = {obj.Neurons.preferredStimulus};
case 'maxRate'
    varargout = {obj.Neurons.maxRate};
case 'backgroundRate'
    varargout = {obj.Neurons.backgroundRate};
case 'integrationTime'
	varargout = {obj.Neurons.integrationTime};
case 'variability'
	varargout = {obj.Neurons.variability};
case 'distribution'
	varargout = {obj.Neurons.distribution};
case 'alpha'
	varargout = {obj.Neurons.alpha};
case 'a'
	varargout = {obj.Neurons.a};
case 'add'
	varargout = {obj.Neurons.add};
case 'R'
	varargout = {obj.Neurons.R};
case 'width'
	varargout = {obj.width};
otherwise
    error([prop_name,' is not a valid Neurons property'])
end
