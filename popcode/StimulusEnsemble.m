classdef StimulusEnsemble
    
    properties
        circular = false;
        ensemble = [];
        width = [];
        pS = [];
        lowerLimit = 0;
        upperLimit = 0;
    end
    
    methods		
        function obj = StimulusEnsemble(varargin)
            switch nargin
            case 3
				switch varargin{1}
                case 'circular'
                    modulo = double(varargin{2});
                    number = double(varargin{3});
                    spacing = modulo / number;
                    obj.circular = modulo;
                    obj.ensemble = [-modulo/2 + spacing : spacing : modulo/2];
                    obj.width = spacing * ones(1, number);
                    obj.pS = 1.0 ./ double(obj.n) .* ones(1, obj.n);
                    
                    obj.lowerLimit = -modulo/2;
                    obj.upperLimit = modulo/2;
				otherwise
					error([varargin{1} ' is not a valid option with three args'])
				end

            case 4
                bottom = double(varargin{2});
                top = double(varargin{3});
                number = floor(varargin{4});

                switch varargin{1}
                case 'linear'
                    obj.circular = false;
                    spacing = (top - bottom) / double(number - 1);
                    obj.ensemble = [bottom : spacing : top];
                    obj.width = diff(obj.ensemble);
                    obj.width = 0.5 * ([obj.width(1) obj.width] + [obj.width obj.width(end)]);
                    
                    obj.lowerLimit = bottom;
                    obj.upperLimit = top;

                otherwise
                    error([strvarargin{1} ' is not a valid option'])
                end

                obj.pS = 1.0 ./ double(obj.n) .* ones(1, obj.n);

            otherwise
                error('Wrong number of arguments')
            end
		end
		
		
		function d = dimensionality(obj)
			d = size(obj.ensemble, 1);
		end
		
		
		function num = n(obj)
			num = size(obj.ensemble, 2);
		end
		
		
		function retStr = char(obj)
			circularStr = num2str(uint32(obj.circular));
			ensembleStr = num2str(obj.ensemble);
			widthStr = num2str(obj.width);
			psStr = num2str(obj.pS);
			
			formatStr = ['Dimensionality:  %.0f\n'...
						 'Modulo:          ' circularStr '\n'...
						 'Stimulus values: ' ensembleStr '\n'...
						 'Bin widths:      ' widthStr '\n'...
						 'P(s):            ' psStr];

			retStr = sprintf(formatStr, obj.dimensionality);
		end
		
	
		function h = entropy(obj)
			h = -sum(obj.pS .* log2(obj.pS));
        end
        
        
        function p = pSint(obj, s)
            % piecewise linear interpolation
            pS = [obj.pS(end) obj.pS];
            ens = [obj.ensemble(1) - 1, obj.ensemble];
            p = interp1q(ens(:), pS', s(:))';
        end
        
	end
end
