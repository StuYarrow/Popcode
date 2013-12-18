classdef OnlineStats < handle
    % OnlineStats Online computation of mean, variance and SEM
    %   
    % Implements the algorithm described in:
    % Donald E. Knuth (1998). The Art of Computer Programming, volume 2: Seminumerical Algorithms, 3rd edn., p. 232. Boston: Addison-Wesley.
    
    properties
        samples = 0;
        runMean = 0;
        runM2 = 0;
        iter = 0;
        logging = true;
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
                    sz = [1 1];
                    maxiter = 1;
                case 2
                    sz = varargin{1};
                    maxiter = varargin{2};
                case 3
                    sz = varargin{1};
                    maxiter = varargin{2};
                    logg = varargin{3};
                    
                    if islogical(logg) && isscalar(logg)
                        obj.logging = logg;
                    else
                        error('Logging option must be scalar boolean')
                    end
                otherwise
                    error('Wrong number of arguments')
            end
            
            assert(length(sz) == 2, 'size must be 2-element row vector')
            
            obj.samples = zeros([sz maxiter]);
            obj.runMean = zeros(sz);
            obj.runM2 = zeros(sz);
            obj.iter = 0;
        end
        
        
        function appendSample(obj, sample)
            % Increment iteration counter
            obj.iter = obj.iter + 1;
            
            % Store sample
            if obj.logging
                obj.samples(:,:,obj.iter) = sample;
            end
            
            % Update mean
            delta = sample - obj.runMean;
            obj.runMean = obj.runMean + delta ./ obj.iter;
            
            % Update second moment and variance
            obj.runM2 = obj.runM2 + delta .* (sample - obj.runMean);
        end
        
        function samp = sample(obj, i)
            samp = obj.samples(:,:,i);
        end
        
        function trim(obj)
            obj.samples = obj.samples(:,:,1:obj.iter);
        end
        
        function mu = mean(obj)
            if obj.logging
                mu = mean(obj.samples(:,:,1:obj.iter), 3);
            else
                mu = obj.runMean;
            end
        end
        
        function se = runSEM(obj)     
            if obj.iter > 1
                se = sqrt((obj.runM2 ./ (obj.iter - 1)) ./ obj.iter);
            else
                se = zeros(size(obj.runMean));
            end
        end
        
        function se = sem(obj)
            if obj.logging
                se = sqrt(var(obj.samples(:,:,1:obj.iter), 0, 3) ./ obj.iter);
            else
                se = obj.runSEM;
            end
        end
        
        function delta = runDelta(obj)
            delta = abs(obj.runSEM ./ obj.runMean);
        end
        
        function objR = minus(objA, objB)
            assert(isa(objA, 'OnlineStats') && isa(objB, 'OnlineStats'), 'Only subtraction of two OnlineStats objects is possible')
            assert(all(size(objA.samples) == size(objB.samples)), 'Instances to be subtracted must have equal dimensionality and iteration count')
            assert(objA.iter == objB.iter, 'Instances to be subtracted must have equal iteration counts')
            assert(objA.logging && objB.logging, 'Instances to be subtracted must have stored samples')
            
            objR = OnlineStats();
            objR.samples = objA.samples - objB.samples;
            objR.iter = objA.iter;
        end
        
        function objR = plus(objA, objB)
            assert(isa(objA, 'OnlineStats') && isa(objB, 'OnlineStats'), 'Only addition of two OnlineStats objects is possible')
            assert(all(size(objA.samples) == size(objB.samples)), 'Instances to be added must have equal dimensionality and iteration count')
            assert(objA.iter == objB.iter, 'Instances to be added must have equal iteration counts')
            assert(objA.logging && objB.logging, 'Instances to be added must have stored samples')
            
            objR = OnlineStats();
            objR.samples = objA.samples + objB.samples;
            objR.iter = objA.iter;
        end
        
    end
    
end

