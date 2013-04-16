classdef LinDiscStimulusEnsemble < DiscreteStimulusEnsemble
    
    properties
        circular = false;
    end
    
    methods
        
        function obj = LinDiscStimulusEnsemble(varargin)
            % superclass constructor
            obj = obj@DiscreteStimulusEnsemble();
            
            switch nargin
            case 0
                % do nothing
                
            case 4
                bottom = double(varargin{2});
                top = double(varargin{3});
                interval = double(varargin{4});
                shift = double(varargin{5});

                assert(top > bottom, 'Upper limit must be greater than lower limit')
                assert(interval > 0 && interval <= (top - bottom), 'Invalid interval value')
                assert(mod(top - bottom, interval) == 0, 'Interval must be a divisor of the range')
                assert(shift >= 0 && shift < interval, 'Invalid shift value')
                
                obj.ensemble = (bottom + shift) : interval : top;
                obj.pS = 1.0 ./ double(obj.n) .* ones(1, obj.n);

            otherwise
                error('Wrong number of arguments')
            end                

        end
        
	end
end
