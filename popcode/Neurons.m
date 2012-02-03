classdef Neurons
	
	properties
		dimensionality = 1;
		popSize = 1;
		preferredStimulus = 0;
		integrationTime = 0;
		distribution = 'Gaussian';
		a = 1.0;
		alpha = 0.5;
		R = [];
		add = 0.0;
		exponent = 1.0;
		truncate = false;
	end
	
	methods
		
		function obj = Neurons(varargin)
		%	NEURONS Neuron population class constructor.
		%	n = Neurons(dimensionality, preferredStimuli, integrationTime, variabilityScheme, variabilityOpts)
		%
		%	dimensionality - stimulus dimensionality (only 1-D stimuli currently supported)
        %   preferred Stimuli - column vector of (N = popsize) characteristic stimuli
		%	integrationTime - spike counting time per trial
		%	variabilityScheme - the type of variablilty model
        %   variabiltyOpts - an array of arguments specific to the chosen variability model
        %
            
			switch nargin
			case 5
				% Standard constructor	
				if length(varargin{1}) == 1 && isnumeric(varargin{1})
					obj.dimensionality = varargin{1};
				else
					error([inputname(1) ' is not a valid stimulus dimensionality'])
				end
				
				if isnumeric(varargin{2}) && size(varargin{2}, 1) == obj.dimensionality
					obj.preferredStimulus = varargin{2}';
					obj.popSize = size(varargin{2}, 2);
				else
					error([inputname(2) ' is not a valid preferred stimulus value or vector'])
				end

				%if length(varargin{3}) == 1 && isnumeric(varargin{3})
				%	obj.maxRate = double(varargin{3}(ones(obj.popSize, 1)));
				%elseif length(varargin{3}) == obj.popSize && isvector(varargin{3}) && isnumeric(varargin{3})
				%	obj.maxRate = reshape(double(varargin{3}), obj.popSize, 1);
				%else
				%	error([inputname(3) ' is not a valid maximum firing rate value or vector for population size ' obj.popSize])
				%end

				%if length(varargin{4}) == 1 && isnumeric(varargin{4})
				%	obj.backgroundRate = double(varargin{4}(ones(obj.popSize, 1)));
				%elseif length(varargin{4}) == obj.popSize && isvector(varargin{4}) && isnumeric(varargin{4})
				%	obj.backgroundRate = reshape(double(varargin{4}), obj.popSize, 1);
				%else
				%	error([inputname(4) ' is not a valid background firing rate value or vector for population size ' obj.popSize])
				%end

				if length(varargin{3}) == 1 && isnumeric(varargin{3})
					obj.integrationTime = double(varargin{3});
				else
					error([inputname(3) ' is not a valid integration time'])
				end

				switch lower(varargin{4})
                    case 'poisson'
                        obj.distribution = 'Poisson';
                        obj.a = [];
                        obj.alpha = [];
                        obj.R = [];
                        obj.add = 0.0;
                        
                    case 'gaussian-independent'
                        obj.distribution = 'Gaussian';
                        obj.a = varargin{5}(1); % need checks
                        obj.alpha = varargin{5}(2); % need checks
                        obj.R = eye(obj.popSize);

                        if length(varargin{5}) == 3
                            obj.add = varargin{5}(3); % need checks
                        else
                            obj.add = 0.0;
                        end
                        
                    case 'gaussian-uniform'
                        obj.distribution = 'Gaussian';
                        obj.a = varargin{5}(1); % need checks
                        obj.alpha = varargin{5}(2); % need checks
                        obj.R = varargin{5}(3) * ~eye(obj.popSize) + eye(obj.popSize);
                        obj.add = 0.0;
                        
                    case 'gaussian-exponential'
                        obj.distribution = 'Gaussian';
                        obj.a = varargin{5}(1); % need checks
                        obj.alpha = varargin{5}(2); % need checks
                        c = varargin{5}(3); % need checks
                        rho = varargin{5}(4); % need checks
                        prefDiff = repmat(obj.preferredStimulus, 1, obj.popSize);
                        prefDiff = prefDiff - prefDiff.';
                        obj.R = c .* exp(-abs(double(prefDiff)) ./ rho) .* ~eye(obj.popSize) + eye(obj.popSize);
                        obj.add = 0.0;
                        
                    case 'gaussian-gaussian'
                        obj.distribution = 'Gaussian';
                        obj.a = varargin{5}(1); % need checks
                        obj.alpha = varargin{5}(2); % need checks
                        c = varargin{5}(3); % need checks
                        beta = 1.0 ./ degToRad(varargin{5}(4)).^2; % need checks
                        prefDiff = repmat(obj.preferredStimulus, 1, obj.popSize);
                        prefDiff = prefDiff - prefDiff.';
                        obj.R = c .* exp((cosd(double(prefDiff)) - 1) .* beta) .* ~eye(obj.popSize) + eye(obj.popSize);
                        obj.add = 0.0;
                        
                    case 'cercal'
                        obj.distribution = 'Gaussian';
                        obj.add = varargin{5}(1);
                        obj.a = varargin{5}(2);
                        obj.alpha = 0.5;
                        obj.R = eye(obj.popSize);
                        obj.exponent = 2.0;
                        
                    otherwise
                        error([varargin{4} ' is not a valid variability regime'])
				end

			otherwise
				error('Wrong number of arguments')
			end
		end
		
		function varargout = mi(obj, method, stim, tol, maxiter)

			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(3) ' is not a SimulusEnsemble object'])
			end

			% obj.popSize x stim.n
			rMean = obj.integrationTime .* meanR(obj, stim);
			rMeanCell = squeeze(mat2cell(rMean, obj.popSize, ones(stim.n, 1)));
            
            % Do distribution-specific one-time prep
            switch obj.distribution
                case 'Gaussian'
                    % Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
                    QCell1 = obj.Q(rMeanCell);

                    % Compute lower triangular Chol(Q) for sampling
                    cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);

                    % Compute upper triangular Chol(Q^-1) for fast PDF computation
                    cholInvQCell = cellfun(@(q) chol(inv(q)), QCell1, 'UniformOutput', false);
                    clear QCell1

                    % Define function for multivariate gaussian sampling
                    % Multiply by Cholesky decomposition of cov matrix Q, and add in mean
                    if obj.truncate
                        fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
                    else
                        fRand = @(m, c, z) m + c * z; % don't truncate
                    end
                    
                case 'Poisson'
                    % Nothing to do here
                    
                otherwise
                    error('Unsupported distribution: %s', obj.distribution)
            end
            
			iter = 0; % iteration counter
            miEst = OnlineStats(1, maxiter);
            adaptive = false;
            tic
            
            cpS = cumsum(stim.pS);
            
            cont = true;
			while cont
                iter = iter + 1;
                
                % Display progress every 100 iterations
				if ~mod(iter, 100)
					fprintf('mi()  iter: %d  val: %.4g  rel. error: %.4g\n', iter, miEst.runMean, miEst.runDelta)
				end

				switch method
				case 'randMC'
					% Sample s from stimulus distribution
					[dummy, bin] = histc(rand(), cpS);
					bin = bin + 1;
					%s = double(stim.ensemble);
					%s = s(bin);

					% Sample r from response distribution
                    switch obj.distribution
                        case 'Gaussian'
                            % Generate vector of independent normal random numbers (mu=0, sigma=1)
                            z = randn(obj.popSize, 1);
                            % Multiply by Cholesky decomposition of cov matrix Q, and add in mean
                            % !!! NOTE NEGATIVE RESPONSES MAY BE TRUNCATED TO ZERO depending on value of obj.truncate !!!
                            r = fRand(rMeanCell{bin}, cholQ{bin}, z);
                        case 'Poisson'
                            % Sample from Poisson distributions
                            r = poissrnd(rMeanCell{bin});
                    end
                    otherwise
                        error('Unsupported method: %s', method)
                end
                
                if ~adaptive
                    % log P(r|s)
                    % Replicate to form a stim.n x stim.n cell array of response vectors
                    rCell = repmat({r}, [stim.n 1]);
                    % Calculate response log probability densities
                    switch obj.distribution
                        case 'Gaussian'
                            lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCell, rMeanCell', cholInvQCell', {'inv'}));
                        case 'Poisson'
                            lpRgS = cell2mat(cellsxfun(@(x, l) sum(poisspdfln(x, l)), rCell, rMeanCell'));
                    end

                    % log P(r,s')
                    % Mutiply P(r|s) and P(s) to find joint distribution
                    pS = stim.pS;
                    lpRS = lpRgS + log(pS');

                    % log P(r)
                    % Calculate marginal by summing over s'
                
                    lpR = logsumexp(lpRS);
                    lpR_sparse = mean([logsumexp(lpRS(1:2:end) + log(2)), logsumexp(lpRS(2:2:end) + log(2))]);
                    
                    % log p(r,s)
                    lpRS = lpRS(bin);

                    if abs((lpR - lpR_sparse) / lpR) > tol
                        % One-shot trapezoid rule is insufficiently accurate; switch to adaptive method
                        adaptive = true;
                    end
                end
                                
                if adaptive
                    pR = [0 0 0];
                    trace = false; % debug flag
                    fAdInt = @quad; % Use the quad function
                    
                    % log p(r,s)
                    lpRS = log(obj.fpSR(stim.ensemble(bin), r, stim));
                    
                    quadTol = exp(lpRS) * tol;
                    
                    switch bin
                        case 1 % Bottom bin - do first 2 bins, remainder
                            [pR(1), fcnt(1)] = fAdInt(@(s) obj.fpSR(s, r, stim), stim.lowerLimit, stim.ensemble(bin+1), quadTol, trace);
                            [pR(2), fcnt(2)] = fAdInt(@(s) obj.fpSR(s, r, stim), stim.ensemble(bin+1), stim.upperLimit, quadTol, trace);
                        case stim.n % Top bin - do first bin, last bin, remainder
                            [pR(1), fcnt(1)] = fAdInt(@(s) obj.fpSR(s, r, stim), stim.lowerLimit, stim.ensemble(1), quadTol, trace);
                            [pR(2), fcnt(2)] = fAdInt(@(s) obj.fpSR(s, r, stim), stim.ensemble(1), stim.ensemble(end-1), quadTol, trace);
                            [pR(3), fcnt(3)] = fAdInt(@(s) obj.fpSR(s, r, stim), stim.ensemble(end-1), stim.upperLimit, quadTol, trace);
                        otherwise % Other bins - do one bin either side, remainder above, remainder below
                            [pR(1), fcnt(1)] = fAdInt(@(s) obj.fpSR(s, r, stim), stim.lowerLimit, stim.ensemble(bin-1), quadTol, trace);
                            [pR(2), fcnt(2)] = fAdInt(@(s) obj.fpSR(s, r, stim), stim.ensemble(bin-1), stim.ensemble(bin+1), quadTol, trace);
                            [pR(3), fcnt(3)] = fAdInt(@(s) obj.fpSR(s, r, stim), stim.ensemble(bin+1), stim.upperLimit, quadTol, trace);
                    end
                    
                    lpR = log(sum(pR));
                    
                    
                end
                
				% log P(s)
				lpS = log(pS(bin));
                                
                % sample MI in bits (convert from log_e to log_2)
                miEst.appendSample((lpRS - (lpR + lpS)) ./ log(2));
                
                % Test halting criteria (SEM, max iterations limit)
				cont = miEst.runDelta > tol & iter < maxiter;
                
                % Impose minimum iteration limit so we get a sensible estimate of SEM
                 cont = cont | iter < 1000;
            end
            
            % Trim unused samples from buffer
            miEst.trim;
            
            % Recompute MI, SEM cleanly
            i = miEst.mean;
            iSem = miEst.sem;
            
            % Report final values
            fprintf('mi() halting  iter: %d  val: %.4g  SEM: %.4g\n', iter, i, iSem)
            
            switch nargout
            case 1
                varargout = {i};
            case 2
                varargout = {i iSem};
            case 3
                varargout = {i iSem miEst};
            otherwise
                error('Unsupported number of return values.')
            end
		end
		
		function varargout = ssiss(obj, n, method, stim, stimOrds, tol, maxiter, timeout)
            
			try
				% Test sanity of neuron indices
				obj.preferredStimulus(n);
			catch err
				error([inputname(2) ' is not a valid neuron index'])
			end
			
			try
				% Test sanity of stimulus ordinate indices
				stim.ensemble(stimOrds);
			catch err
				error([inputname(5) ' is not a valid stimulus ordinate index'])
			end
			
			if ~any(strcmp(method, {'quadrature' 'randMC' 'quasirandMC'}))
				error([method ' is not a valid SSI calculation method'])
            end

            if ~isa(stim, 'StimulusEnsemble')
				error([inputname(4) ' is not a SimulusEnsemble object'])
            end
			
			% Create mask for calculating specific stimulus ordinates only
			if ~isempty(stimOrds)
				sMask = false(stim.n, 1);
				sMask(stimOrds) = true;
				sMaskN = sum(sMask + 0);
			else
				sMask = true(stim.n, 1);
				sMaskN = stim.n;
				stimOrds = 1:stim.n;
			end
			
			% Get mean responses for each stimulus ordinate
			% obj.popSize x stim.n
			rMean = obj.integrationTime .* obj.meanR(stim);
			rMeanCell = squeeze(mat2cell(rMean, obj.popSize, ones(stim.n, 1)));
            
            % Setup for computing marginals
            if ~isempty(n)
                % Create logical vector (mask) identifying neurons that are *not* part of the marginal SSI
                margMask = true(obj.popSize, 1);
                margMask(n) = false;

                % Get mean responses for each stimulus ordinate
                rMeanMargCell = cellfun(@(r) r(margMask), rMeanCell, 'UniformOutput', false);
            end
            
            switch obj.distribution
                case 'Gaussian'
                    % Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
                    QCell1 = obj.Q(rMeanCell);

                    % Compute lower triangular Chol(Q) for sampling
                    cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);

                    % Compute upper triangular Chol(Q^-1) for fast PDF computation
                    cholInvQCell = cellfun(@(q) chol(inv(q)), QCell1, 'UniformOutput', false);
                    
                    if ~isempty(n)
                        % Compute mean response dependent cov matrix stack Q
                        QCellMarg1 = cellfun(@(q) q(margMask, margMask), QCell1, 'UniformOutput', false);

                        % Invert Q matrices and compute Cholesky decomps
                        cholInvQCellMarg = cellfun(@(q) chol(inv(q)), QCellMarg1, 'UniformOutput', false);
                        clear QCellMarg1
                    end
                    
                    clear QCell1
                    
                    % Define function for multivariate gaussian sampling
                    % Multiply by Cholesky decomposition of cov matrix Q, and add in mean
                    if obj.truncate
                        fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
                    else
                        fRand = @(m, c, z) m + c * z; % don't truncate
                    end
                    
                case 'Poisson'
                    % Nothing to be done here
                    
                otherwise
                    error('Unsupported distribution: %s', obj.distribution)
            end
			
            % Initialise main loop, preallocate MC sample arrays
			iter = 0;
			cont = true;
            
            Issi = OnlineStats(sMaskN, maxiter);
            Isur = OnlineStats(sMaskN, maxiter);
            IssiMarg = OnlineStats(sMaskN, maxiter);
            IsurMarg = OnlineStats(sMaskN, maxiter);
			
            % Main MC sampling loop
            while cont
                iter = iter + 1;

				if ~mod(iter, 10)
					fprintf('SSISS iter: %d of %d, rel. error: %.4g\n', iter, maxiter, mean(Issi.runDelta))
				end

				switch method
                    case 'randMC'
                        % Sample r from response distribution
                        switch obj.distribution
                            case 'Gaussian'
                                % Generate vector of independent normal random numbers (mu=0, sigma=1)
                                zCell = mat2cell(randn(obj.popSize, sMaskN), obj.popSize, ones(sMaskN, 1));
                                % Multiply by Cholesky decomposition of cov matrix Q, and add in mean
                                % !!! NOTE NEGATIVE RESPONSES MAY BE TRUNCATED TO ZERO, SEE ABOVE !!!
                                rCell = cellfun(fRand, rMeanCell(sMask), cholQ(sMask), zCell, 'UniformOutput', false); % stim.n cell array of obj.popSize vectors
                                
                            case 'Poisson'
                                % Sample from Poisson distributions
                                rCell = cellfun(@poissrnd, rMeanCell(sMask), 'UniformOutput', false);
                        end
                        
                    otherwise
                        error('Unsupported method: %s', method)
				end

				% log P(r|s')
				% Calculate response probability densities
                switch obj.distribution
                    case 'Gaussian'
                        lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCell, rMeanCell', cholInvQCell', {'inv'}));
                    case 'Poisson'
                        lpRgS = cell2mat(cellsxfun(@(x, l) sum(poisspdfln(x, l)), rCell, rMeanCell'));
                end
                
				% log P(r,s')
				% Mutiply P(r|s) and P(s) to find joint distribution
				lpRS = bsxfun(@plus, lpRgS, log(stim.pS')); % stim'.n x stim

				% log P(r)
				% Calculate marginal by summing over s'
				lpR = logsumexp(lpRS, 1);
				
				% log P(s'|r)
				% Divide joint by marginal P(r)
				lpSgR = bsxfun(@minus, lpRS, lpR);

				% H(s'|r), in bits, converting from log_e to log_2
				hSgR = -sum(exp(lpSgR) .* (lpSgR ./ log(2)), 1);

				% Sample specific information Isp(r)
				% Specific information; reduction in stimulus entropy due to observation of r
				Issi.appendSample(stim.entropy - hSgR);
                
                % Sample specific surprise
				% log_2( P(r|s) / P(r) )
				% Accumulate samples
                Isur.appendSample((diag(lpRgS(stimOrds,:))' - lpR) ./ log(2));
				
                if exist('margMask', 'var')
					% If we are calculating a marginal SSI, compute the SSI for remaining neurons
                    
					% Mask out neurons of interest in response vectors
					rCellMarg = cellfun(@(r) r(margMask), rCell, 'UniformOutput', false);
                    
					% log P(r|s)
                    switch obj.distribution
                        case 'Gaussian'
                            lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCellMarg, rMeanMargCell', cholInvQCellMarg', {'inv'}));
                        case 'Poisson'
                            lpRgS = cell2mat(cellsxfun(@(x, l) sum(poisspdfln(x, l)), rCellMarg, rMeanMargCell'));
                    end
                    
					% log P(r,s)
					% Multiply P(r|s) and P(s) to find joint distribution
					lpRS = bsxfun(@plus, lpRgS, log(stim.pS')); % stim'.n x stim
                    
					% log P(r)
					% Calculate marginal by summing over s'
					lpR = logsumexp(lpRS, 1);
                    
					% log P(s|r)
					% Divide joint by marginal P(r)
					lpSgR = bsxfun(@minus, lpRS, lpR);
                    
					% H(s|r), in bits, converting from log_e to log_2
					hSgR = -sum(exp(lpSgR) .* (lpSgR ./ log(2)), 1);
                    
					% Isp(r)
					% Specific information; reduction in stimulus entropy due to observation of r
                    % Compute MC sample
                    IssiMarg.appendSample(stim.entropy - hSgR);
					
                    % Specific surprise
					% log_2( P(r|s) / P(r) )
					% Compute MC sample
                    IsurMarg.appendSample((diag(lpRgS(stimOrds,:))' - lpR) ./ log(2));
                end
                
                % Test halting criteria (SEM, max iterations limit, timeout)
                if exist('margMask', 'var')
                    delta = mean([Issi.runDelta IssiMarg.runDelta]);
                else
                    delta = mean(Issi.runDelta);
                end
                
                cont = delta > tol & iter < maxiter;
                
                % If the wall clock is running, check the elapsed time
                try
                    cont = cont & toc < timeout;
                catch
                end
                
                % Impose minimum iteration limit so we get a valid estimate of SEM
                cont = cont | iter < 10;
            end
            
            try
                fprintf('SSISS iter: %d  elapsed time: %.4f seconds\n', iter, toc)
            catch
                fprintf('SSISS iter: %d\n', iter)
            end
            
            % Trim sample arrays
            Issi.trim;
            Isur.trim;
            IssiMarg.trim;
            IsurMarg.trim;
                        
            % Recalculate means cleanly
            fullSSI = Issi.mean;
            fullIsur = Isur.mean;
            
            if exist('margMask', 'var')
                remSSI = IssiMarg.mean;
                SSI = fullSSI - remSSI;
                remIsur = IsurMarg.mean;
            else
                remSSI = fullSSI;
                SSI = fullSSI;
                remIsur = fullIsur;
            end
            
			switch nargout
			case 1
				varargout = {SSI};
			case 2
				varargout = {SSI iter};
			case 3
				varargout = {fullSSI remSSI iter};
            case 4
                varargout = {fullSSI remSSI iter Issi IssiMarg};
			case 5
				varargout = {fullSSI remSSI fullIsur remIsur iter};
            case 9
                varargout = {fullSSI remSSI fullIsur remIsur iter Issi Isur IssiMarg IsurMarg};
			otherwise
				error('Unsupported number of outputs')
			end
        end
		
		function varargout = fisher(obj, method, stim, tol, maxiter)
            % Wrapper function for calculating Fisher information
            
            switch method
                case 'analytic'
                    switch obj.distribution
                        case 'Gaussian'
                            [J_mean, J_cov] = obj.fisher_analytic_gauss(stim);
                            
                            switch nargout
                                case 0
                                    varargout = {J_mean + J_cov};
                                case 1
                                    varargout = {J_mean + J_cov};
                                case 2
                                    varargout = {J_mean J_cov};
                                otherwise
                                    error('Wrong number of outputs')
                            end
                            
                        case 'Poisson'
                            J = obj.fisher_analytic_poiss(stim);
                            
                            if nargout == 1
                                varargout = {J};
                            else
                                error('Wrong number of outputs')
                            end
                            
                        otherwise
                            error('Unsupported distribution: %s', obj.distribution)
                    end
                    
                case 'randMC'
                    [J samples] = obj.fisher_mc(method, stim, tol, maxiter);
                    
                    switch nargout
                        case 1
                            varargout = {J};
                        case 2
                            varargout = {J samples};
                        otherwise
                            error('Wrong number of outputs')
                    end
                    
                otherwise
                    error('"%s" is not a valid FI calculation method', method)
            end
		end
		
		function [fullFisher remainderFisher] = margFisher(obj, nMarg, method, stim, tol)
            % Function for calculating fisher of population with and
            % without neuron(s) of interest
            
			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(4) ' is not a SimulusEnsemble object'])
			end

			% Compute FI of full population
			fullFisher = fisher(obj, method, stim, tol);
			
			% Create population of remaining neurons
			remainingNeurons = obj.remove(nMarg);
			
			% FI of remaining neurons
			remainderFisher = fisher(remainingNeurons, method, stim, tol);
		end
		
		function ifish = Ifisher(obj, stim)
            % Function for computing I_Fisher, see:
            %
            % Brunel N, Nadal J (1998)
            % Mutual information, Fisher information, and population coding.
            % Neural Comput 10:1731?1757.
            
            if false
                % Fisher information
                fish = obj.fisher('analytic', stim, 0);

                % Ignore zero Fisher information values
                zeroJ = find(fish == 0);
                pS = stim.pS;

                if ~isempty(zeroJ)
                    fish(zeroJ) = [];
                    pS(zeroJ) = [];
                    pS = pS ./ sum(pS);
                end
                
                ifish = stim.entropy - 0.5 .* sum(pS .* (1.0 + log2(pi) + log2(exp(1)) - log2(fish)));
            else
                ifish = stim.entropy - 0.5 .* quad(@(s) stim.pSint(s) .* (1.0 + log2(pi) + log2(exp(1)) - log2(obj.fisher('analytic', s, 0))), stim.ensemble(1) - diff(stim.ensemble(1:2)), stim.ensemble(end));
            end
		end
		
		function varargout = mIfisher(obj, nMarg, stim)
            % Function for computing the marginal I_Fisher
            
            % Get whole-population I_Fisher
            ifishFull = obj.Ifisher(stim);
            
            % Create population of remaining neurons
			remainingNeurons = obj.remove(nMarg);
            
            % Compute remainder I_Fisher
            ifishRem = remainingNeurons.Ifisher(stim);
            
            switch nargout
                case 1
                    varargout = {ifishFull - ifishRem};
                case 2
                    varargout = {ifishFull ifishRem};
                otherwise
                    error('Wrong number of outputs')
            end
		end
		
		function ssif = SSIfisher(obj, n, fisherMethod, stim, tol)
            % Function for computing SSI_Fisher, see:
            %
            % Yarrow S, Challis E, Series P (2012)
            % Fisher and Shannon information in finite neural populations.
            % Neural Comput (in review).
		
			% S stimulus
			sMat = repmat(stim.ensemble, [stim.n 1]);
			% sHat stimulus estimate
			sHatMat = repmat(stim.ensemble', [1 stim.n]);
			% log p(S)
			psMat = repmat(stim.pS, [stim.n 1]);
			
			if stim.circular
				% circular difference
				dS = mod(sHatMat - sMat, stim.circular);
				i = find(dS > (0.5 * stim.circular));
				dS(i) = dS(i) - stim.circular;
				i = find(dS < (-0.5 * stim.circular));
				dS(i) = dS(i) + stim.circular;
			else
				% linear difference
				dS = sHatMat - sMat;
			end
			
			if ~isempty(n)
				[fullFI remFI] = obj.margFisher(n, fisherMethod, stim, tol);
				
				% Compute SSIfisher excluding cells of interest
				% sigma(s)
				% Compute SD of optimal estimator as a function of the stimulus
				sigma = remFI .^ -0.5;
				sigmaMat = repmat(sigma, [stim.n 1]);

				% log p(sHat|S)
				lpsHat_s = cellfun(@mvnormpdfln, num2cell(dS), num2cell(zeros([stim.n stim.n])), num2cell(sigmaMat));
				
				% log p(S)
				lpS = log(psMat);
				% log p(sHat,S)
				lpsHats = lpsHat_s + lpS;
				% log p(sHat)
				lpsHat = logsumexp(lpsHats, 2);
				% log p(S|sHat)
				lps_sHat = lpsHats - repmat(lpsHat, [1 stim.n]);

				% H(s|sHat) as a function of sHat
				Hs_sHat = -sum(exp(lps_sHat) .* lps_sHat ./ log(2), 2);

				% Isp(sHat) specific information
				isp = stim.entropy - Hs_sHat;

				% SSIfisher
                ssifRem = sum(exp(lpsHat_s + repmat(log(isp), [1 stim.n])), 1);
			else
				fullFI = obj.fisher(fisherMethod, stim, tol);
			end
			
			% Compute SSIfisher for full population
			% sigma(s)
			% Compute SD of optimal estimator as a function of the stimulus
			sigma = fullFI .^ -0.5;
			sigmaMat = repmat(sigma, [stim.n 1]);

			% log p(sHat|S)
			lpsHat_s = cellfun(@mvnormpdfln, num2cell(dS), repmat({0}, [stim.n stim.n]), num2cell(sigmaMat));
			
			% log p(S)
			lpS = log(psMat);
			% log p(sHat,S)
			lpsHats = lpsHat_s + lpS;
			% log p(sHat)
			lpsHat = logsumexp(lpsHats, 2);
			% log p(S|sHat)
			lps_sHat = lpsHats - repmat(lpsHat, [1 stim.n]);

			% H(s|sHat) as a function of sHat
			Hs_sHat = -sum(exp(lps_sHat) .* lps_sHat ./ log(2), 2);

			% Isp(sHat) specific information
			isp = stim.entropy - Hs_sHat;
			
			% SSIfisher
            ssifFull = sum(exp(lpsHat_s + repmat(log(isp), [1 stim.n])), 1);
			
			if ~isempty(n)
				ssif = ssifFull - ssifRem;
			else
				ssif = ssifFull;
			end
		end
				
		function q = Q(obj, resp)
			q = cellfun(@(r) (obj.add .* obj.R + obj.a .* obj.R .* (r * r').^obj.alpha).^obj.exponent, resp, 'UniformOutput', false);
		end
		
		function retStr = char(obj)
			%	CHAR Text representation of Neurons object
			%	aString = char(object)

			if ~isscalar(obj)
				retStr = 'Array of Neurons objects';
				return
			end

			prefStr = char(obj.preferredStimulus);
			%maxRateStr = char(obj.maxRate);
			%backgroundRateStr = char(obj.backgroundRate);
			rStr = display(obj.R);

			formatStr = ['Population size: %.0f\n'...
						 'Stimulus dimensionality: %.0f\n'...
						 'Preferred stimuli:\n'...
						  prefStr '\n'...
						 %'Maximum firing rate (Hz):\n'...
						 % maxRateStr '\n'...
						 %'Background firing rate (Hz):\n'...
						 % backgroundRateStr '\n'...
						 'Integration time: %.3f s\n'...
						 'Variability distribution: %s\n'...
						 'Variability coefficient a: %.3f\n'...
						 'Variability exponent alpha: %.3f\n'...
						 'Correlation matrix:\n'...
                         '%s\n'];

			retStr = sprintf(formatStr, obj.popSize, double(obj.dimensionality), obj.integrationTime, obj.distribution, obj.a, obj.alpha, rStr);
		end
		
		function retVal = tcplot(obj, stim)
			r = meanR(obj, stim);
			[stims ind] = sort(double(stim.ensemble));
			retVal = plot(stims, r(:,ind));
		end
		
		function varargout = remove(obj, nMarg)
			if ~isempty(nMarg)
				% Create logical vector (mask) identifying neurons that are *not* part of the marginal SSI
				margMask = ones(obj.popSize, 1);
				margMask(nMarg) = false;
				% Number of remaining neurons
				nMarg = sum(margMask);
				margMask = logical(margMask);
			else
				error('Must specify index of a neuron or neurons')
			end
			
			obj.preferredStimulus = obj.preferredStimulus(margMask);
			obj.popSize = nMarg;
			
			Rmask = logical((margMask+0) * (margMask+0)');
			obj.R = reshape(obj.R(Rmask), [nMarg nMarg]);
			
			switch nargout
			case 1
				varargout = {obj};
			case 2
				varargout = {obj margMask};
			otherwise
				error('Wrong number of outputs')
			end
					
		end
		
    end
    
    
    methods (Access = protected)
        
        function [J_mean, J_cov] = fisher_analytic_gauss(obj, stim)
            % Info=Fisher(N, alpha, correlation_type, param_correlation)
            % compute I_mean, I_cov and I_mean and I_cov if neurons are independent 
            % when the tuning function is circular normal and correlation is toeplitz
            % 
            % PARAMETERS:
            % -----------
            % - N is the number of neurons
            % - f is the tuning curve
            % - f_prime is the derivative of the tuning curve
            % - Q is the covariance matrix
            % - k is the fano factor, now a vector
            % - k prime is its derivative
            %
            % RETURNS:
            % ---------
            % Info(1):  I_mean
            % Info(2):  I_mean_ind if the correlations were set to 0 (same variance)
            % Info(3):  I_cov
            % Info(4):  I_cov_ind if the correlations were set to 0
            % Info(5):  I_mean+I_cov
            % Info(6):  I_SQE

            % pseries@gatsby.ucl.ac.uk  2/20/2005
            % s.yarrow@ed.ac.uk 2008-2011

            f = obj.integrationTime .* obj.meanR(stim);
            f_prime = obj.integrationTime .* obj.dMeanR(stim);

            g0 = f_prime ./ f;

            % ==============================================================
            % Correlation matrix, Covariance matrix, inverse and derivative.
            % Fourier tranform (eigenvalues)
            % ==============================================================
            
            if isa(stim, 'StimulusEnsemble')
                nStim = stim.n;
            else
                nStim = length(stim);
            end
            
            fCell = squeeze(mat2cell(f, obj.popSize, ones(nStim, 1)));
            f_primeCell = squeeze(mat2cell(f_prime, obj.popSize, ones(nStim, 1)));
            g0Cell = squeeze(mat2cell(g0, obj.popSize, ones(nStim, 1)));
            QCell1 = obj.Q(fCell);
            Q_inv = cellfun(@inv, QCell1, 'UniformOutput', false); % inverse
            k = repmat({obj.a(ones(obj.popSize,1))}, [1 nStim]);
            k_prime = repmat({zeros(obj.popSize,1)}, [1, nStim]);

            fQ_prime = @(q, kk, kp, gz) obj.alpha * (diag(gz) * q + q * diag(gz)) + (diag(kp ./ kk)) * q;
            Q_prime = cellfun(fQ_prime, QCell1, k, k_prime, g0Cell, 'UniformOutput', false); % derivative

            % ==============================================================
            % Fisher Information
            % J_mean, J_cov
            % (J = J_mean + J_cov)
            % ==============================================================

            if obj.popSize > 1
                J_mean = cellfun(@(fp, qi) fp' * qi * fp, f_primeCell, Q_inv);         % direct method
                %Info(2) = cellfun(@(fp, di) fp * di * fp', f_primeCell, D_inv);
            else
                J_mean = cellfun(@(fp, q) fp.^2 / q, f_primeCell, QCell1);
            end

            J_cov = cellfun(@(qp, qi) 0.5 * trace(qp * qi * qp * qi), Q_prime, Q_inv);  % direct method for Icov
        end
        
        function J = fisher_analytic_poiss(obj, stim)
            % TC and TC derivative
            f = obj.meanR(stim);
            f_prime = obj.dMeanR(stim);
            
            % Fisher information, see:
            % Bethge M, Rotermund D, Pawelzik K (2002)
            % Optimal short-term population coding: when Fisher information fails.
            % Neural Comput 14:2317?2351.
            %
            J = obj.integrationTime .* sum(f_prime.^2 ./ f, 1);
        end
        
        function [J FI] = fisher_mc(obj, method, stim, tol, maxiter)
            % Generate offset stimuli
            stim2 = stim;
            dTheta = 0.1 .* stim2.width;
            stim2.ensemble = stim2.ensemble + dTheta;

            % obj.popSize x stim.n
            rMean1 = obj.integrationTime .* obj.meanR(stim);
            rMeanCell1 = squeeze(mat2cell(rMean1, obj.popSize, ones(stim.n, 1)));
            rMean2 = obj.integrationTime .* obj.meanR(stim2);
            rMeanCell2 = squeeze(mat2cell(rMean2, obj.popSize, ones(stim2.n, 1)));
            
            % Concat (stim) and (stim + dTheta) mean arrays
            rMeanCellDiff = [rMeanCell1 ; rMeanCell2];
            
            % Do distribution-specific setup
            switch obj.distribution
                case 'Gaussian'
                    % Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
                    QCell1 = obj.Q(rMeanCell1);
                    QCell2 = obj.Q(rMeanCell2);

                    % Compute Cholesky decomps of Q matrices
                    cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);

                    % Invert Q matrices and compute Cholesky decomps
                    cholInvQCell1 = cellfun(@(q) chol(inv(q)), QCell1, 'UniformOutput', false);
                    cholInvQCell2 = cellfun(@(q) chol(inv(q)), QCell2, 'UniformOutput', false);
                    clear QCell1 QCell2

                    % Concat (stim) and (stim + dTheta) chol(Q^-1) arrays
                    cholInvQCellDiff = [cholInvQCell1 ; cholInvQCell2];
                    clear cholInvQCell1 cholInvQCell2

                    % Define function for multivariate gaussian sampling
                    % Multiply by Cholesky decomposition of cov matrix Q, and add in mean
                    if obj.truncate
                        fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
                    else
                        fRand = @(m, c, z) m + c * z; % don't truncate
                    end
                    
                case 'Poisson'
                    % Nothing to do here
                    
                otherwise
                    error('Unsupported distribution: %s', obj.distribution)
            end

            iter = 0;
            cont = true;

            FI = OnlineStats(stim.n, maxiter);

            while cont
                iter = iter + 1;

                if ~mod(iter, 100)
                    fprintf('Fisher iter: %d, rel. error: %.4g\n', iter, mean(FI.runDelta))
                end

                switch method
                    case 'randMC'
                        % Sample r from response distribution
                        switch obj.distribution
                            case 'Gaussian'
                                % Generate vector of independent normal random numbers (mu=0, sigma=1)
                                zCell = mat2cell(randn(obj.popSize, stim.n), obj.popSize, ones(stim.n, 1));
                                % Multiply by Cholesky decomposition of cov matrix Q, and add in mean
                                rCell = cellfun(fRand, rMeanCell1, cholQ, zCell, 'UniformOutput', false); % stim.n cell array of {obj.popSize}
                                
                            case 'Poisson'
                                % Sample from Poisson
                                rCell = cellfun(@poissrnd, rMeanCell1, 'UniformOutput', false);
                        end
                        
                    otherwise
                        error('Unsupported method: %s', method)
                end

                % log P(r|s)
                % Calculate response probability densities
                switch obj.distribution
                    case 'Gaussian'
                        lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCell, rMeanCellDiff, cholInvQCellDiff, {'inv'})); % 1 x stim.n
                    case 'Poisson'
                        lpRgS = cell2mat(cellsxfun(@(x, l) sum(poisspdfln(x, l)), rCell, rMeanCellDiff)); % 1 x stim.n
                end

                % d/ds log P(r|s)
                % Find derivative
                dlpRgS = diff(lpRgS, 1, 1) ./ dTheta;

                % (d/ds log P(r|s))^2
                FI.appendSample(dlpRgS .^ 2);

                % Test halting criteria (SEM, max iterations limit, timeout)
                delta = mean(FI.runDelta);
                cont = all([delta > tol, iter < maxiter]);

                % Impose minimum iteration limit so we get a sensible estimate of SEM
                cont = cont | iter < 10;
            end

            fprintf('FI completed iter: %d\n', iter)

            % Trim sample array
            FI.trim;
            
            % Return output
            J = FI.mean;
        end
        
        function pSR = fpSR(obj, s, r, stim)
            % mean response given s
            rMean = obj.integrationTime .* meanR(obj, s);
            rMeanCell = squeeze(mat2cell(rMean, obj.popSize, ones(length(s), 1)));
            rCell = repmat({r}, [length(s) 1]);
            
            switch obj.distribution
            case 'Gaussian'
                % covariance matrix
                Q = obj.Q(rMeanCell);
                
                % log p(r|s)
                lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCell, rMeanCell', {[]}, Q'))';
            case 'Poisson'
                % log p(r|s)
                lpRgS = sum(poisspdfln(r, rMean));
            end
            
            % log p(s) via linear piecewise interpolation on stimulus distribution
            lpS = log(stim.pSint(s));
            
            % p(s,r) = exp(log p(s) + log p(r|s))
            pSR = exp(lpS + lpRgS);
        end
        
    end

end
