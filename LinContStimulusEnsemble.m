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
                
            case 2
                obj.ensemble = double(varargin{1}(:)');
                obj.pS = double(varargin{2}(:)');
                
                assert(all(diff(diff(obj.ensemble)) < 1e-8), 'Stimulus distribution points must be regularly spaced')
                
                obj.width = diff(obj.ensemble(1:2));
                obj.lowerLimit = obj.ensemble(1);
                obj.upperLimit = obj.ensemble(end);
                obj.pS = obj.pS ./ obj.integrate(obj.pS, 2);
                
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
                obj.pS = 1.0 ./ (top - bottom) .* ones(1, obj.n);
                
            otherwise
                error('Wrong number of arguments')
            end
		end
	
        function p = pSint(obj, s)
            % piecewise linear interpolation
            p = interp1q(obj.ensemble', obj.pS', s(:))';
        end
        
        function integral = integrate(obj, ords, dim)
            % trapezoid rule
            integral = trapz(ords, dim) .* obj.width;
        end
        
	end
end
