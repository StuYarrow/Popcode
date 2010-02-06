function obj = set(obj, varargin)
% SET Set Neurons properties and return the updated object

property_argin = varargin;

while length(property_argin) >= 2,
	prop = property_argin{1};
	val = property_argin{2};
	property_argin = property_argin(3:end);
	
	try
		obj.Neurons = set(obj.Neurons, prop, val);
	catch
		switch prop
		case 'width'
			if (length(val) == 1 | length(val) == obj.Neurons.popSize) & isnumeric(val) & isvector(val)
				obj.width = val;
			else
				error(['width value is not a valid for population size ' num2str(obj.Neurons.popSize)])
			end	
		otherwise
			error('Invalid field name')
		end
	end
end