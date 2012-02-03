classdef OnlineStats < handle
    % OnlineStats Online computation of mean, variance and SEM
    %   
    % Implements the algorithm described in:
    % Donald E. Knuth (1998). The Art of Computer Programming, volume 2: Seminumerical Algorithms, 3rd edn., p. 232. Boston: Addison-Wesley.
    
    properties
        samples = 0;
        runMean = 0;
        runM2 = 0;
        runVar = 0;
        runSEM = 0;
        iter = 0;
    end
    
    methods
        
        function obj = OnlineStats(varargin)
            % OnlineStats class constructor
            % 
            % 
            % Usage obj = OnlineStats(l, maxiter)
            % l is the length of each sample vector
            % maxiter is the max number of iterations, for memory preallocation
            
            switch nargin
                case 0
                    l = 1;
                    maxiter = 1;
                case 2
                    l = varargin{1};
                    maxiter = varargin{2};
                otherwise
                    error('Wrong number of arguments')
            end
            
            obj.samples = zeros(maxiter, l);
            obj.runMean = zeros(1, l);
            obj.runM2 = zeros(1, l);
            obj.iter = 0;
        end
        
        
        function appendSample(obj, sample)
            % Increment iteration counter
            obj.iter = obj.iter + 1;
            
            % Store sample
            obj.samples(obj.iter,:) = sample;
            
            % Update mean
            delta = obj.samples(obj.iter,:) - obj.runMean;
            obj.runMean = obj.runMean + delta ./ obj.iter;
            
            % Update second moment and variance
            obj.runM2 = obj.runM2 + delta .* (obj.samples(obj.iter,:) - obj.runMean);
            obj.runVar = (obj.runM2 ./ (obj.iter - 1));
            
            % Calculate SEM
            obj.runSEM = sqrt(obj.runVar ./ obj.iter);
        end
        
        
        function trim(obj)
            obj.samples = obj.samples(1:obj.iter,:);
        end
        
        
        function mu = mean(obj)
            mu = mean(obj.samples(1:obj.iter,:), 1);
        end
        
        
        function err = sem(obj)
            err = sqrt(var(obj.samples(1:obj.iter,:), 1) ./ obj.iter);
        end
        
        
        function delta = runDelta(obj)
            delta = abs(obj.runSEM ./ obj.runMean);
        end
    end
    
end

