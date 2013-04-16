classdef LinContStimulusEnsemble < ContinuousStimulusEnsemble
    
    properties
        circular = false;
    end
    
    methods
        
        function obj = LinContStimulusEnsemble(varargin)
            % superclass constructor
            obj = obj@ContinuousStimulusEnsemble();
            
            switch nargin
            case 0
                % do nothing
                                
            case 3
                bottom = double(varargin{1});
                top = double(varargin{2});
                number = double(varargin{3});
                spacing = (top - bottom) / double(number - 1);
                
                assert(mod(number, 1) == 0, 'Non-integer number of stimuli')
                assert(top > bottom, 'Upper limit must be greater than lower limit')

                obj.ensemble = bottom : spacing : top;
                obj.width = spacing;
                obj.lowerLimit = bottom;
                obj.upperLimit = top;
                obj.pS = 1.0 ./ double(obj.n) .* ones(1, obj.n);
                
            otherwise
                error('Wrong number of arguments')
            end
		end
	
        function p = pSint(obj, s)
            % piecewise linear interpolation
            p = interp1q(obj.ensemble', obj.pS', s(:))';
        end
        
	end
end
