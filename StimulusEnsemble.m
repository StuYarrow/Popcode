classdef StimulusEnsemble
    
    properties
        ensemble = [];
        pS = [];
    end
    
    methods
        
        function obj = StimulusEnsemble(varargin) 
            assert(nargin == 0, 'Wrong number of arguments to constructor')
        end
        
		
		function d = dimensionality(obj)
			d = size(obj.ensemble, 1);
		end
		
		
		function num = n(obj)
			num = size(obj.ensemble, 2);
        end
        
        
        function h = entropy(obj)
            h = -obj.integrate(obj.pS .* log2(obj.pS), 2);
        end
        
    end
    
    methods (Abstract)
        integrate(obj, dim)
    end
end
