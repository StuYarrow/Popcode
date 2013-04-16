classdef DiscreteStimulusEnsemble < StimulusEnsemble
	
	properties
		continuous = false;
	end
	
	methods

		function obj = DiscreteStimulusEnsemble(varargin)
		% DiscreteStimulusEnsemble/DiscreteStimulusEnsemble Constructor for DiscreteStimulusEnsemble object
		% obj = DiscreteStimulusEnsemble()
        %
            
            assert(nargin == 0, 'Wrong number of arguments')
            
			% Superclass constructor
            obj = obj@StimulusEnsemble();
        end
        
    end
end