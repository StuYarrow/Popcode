classdef CircContStimulusEnsemble < ContinuousStimulusEnsemble
    
    properties
        circular = 1;
    end
    
    methods
        
        function obj = CircContStimulusEnsemble(varargin)
            % superclass constructor
            obj = obj@ContinuousStimulusEnsemble();
            
            switch nargin
            case 0
                % do nothing
                
            case 2
                modulo = double(varargin{1});
                number = double(varargin{2});
                spacing = modulo / number;
                
                assert(mod(number, 1) == 0, 'Non-integer number of stimuli')
                
                obj.circular = modulo;
                obj.ensemble = -modulo/2 + spacing : spacing : modulo/2;
                obj.width = spacing;
                obj.lowerLimit = -modulo/2;
                obj.upperLimit = modulo/2;
                obj.pS = 1.0 ./ modulo .* ones(1, obj.n);
                
            otherwise
                error('Wrong number of arguments')
            end
		end
	
        function p = pSint(obj, s)
            % piecewise linear interpolation
            pS = [obj.pS(end), obj.pS];
            ens = [obj.ensemble(1) - obj.width, obj.ensemble];
            p = interp1q(ens(:), pS', s(:))';
        end
        
        function h = entropy(obj)
			h = -sum(obj.pS .* log2(obj.pS)) .* obj.width;
        end
        
	end
end
