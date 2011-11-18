classdef Neurons
	
	properties
		dimensionality = 1;
		preferredStimulus = 0;
		popSize = 1;
		maxRate = 0;
		backgroundRate = 0;
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
		%	n = Neurons(dimensionality, preferredStimulus, maxRate, backgroundRate, integrationTime, variabilityScheme, variabilityOpts)
		%
		%	dimensionality - stimulus dimensionality (only 1-D stimuli currently supported)
		%	maxRate - maximum firing rate (Hz)
		%	backgroundRate - background (spontaneous) firing rate (Hz)
		%	variability - a Variability object
		%	integrationTime - spike counting time per trial
		%
		%	maxFiringRate and backgroundFiringRate can be scalars or vectors of length popSize

			%varargin
			%size(varargin{2})
                        
			switch nargin
			case 7
				% Standard constructor	
				if length(varargin{1}) == 1 && isnumeric(varargin{1})
					obj.dimensionality = varargin{1};
				else
					error([inputname(1) ' is not a valid stimulus dimensionality'])
				end
				
				if isnumeric(varargin{2}) & size(varargin{2}, 1) == obj.dimensionality
					obj.preferredStimulus = varargin{2}';
					obj.popSize = size(varargin{2}, 2);
				else
					error([inputname(2) ' is not a valid preferred stimulus value or vector'])
				end

				if length(varargin{3}) == 1 && isnumeric(varargin{3})
					obj.maxRate = double(varargin{3}(ones(obj.popSize, 1)));
				elseif length(varargin{3}) == obj.popSize && isvector(varargin{3}) && isnumeric(varargin{3})
					obj.maxRate = reshape(double(varargin{3}), obj.popSize, 1);
				else
					error([inputname(3) ' is not a valid maximum firing rate value or vector for population size ' obj.popSize])
				end

				if length(varargin{4}) == 1 && isnumeric(varargin{4})
					obj.backgroundRate = double(varargin{4}(ones(obj.popSize, 1)));
				elseif length(varargin{4}) == obj.popSize && isvector(varargin{4}) && isnumeric(varargin{4})
					obj.backgroundRate = reshape(double(varargin{4}), obj.popSize, 1);
				else
					error([inputname(4) ' is not a valid background firing rate value or vector for population size ' obj.popSize])
				end

				if length(varargin{5}) == 1 && isnumeric(varargin{5})
					obj.integrationTime = double(varargin{5});
				else
					error([inputname(5) ' is not a valid integration time'])
				end

				switch lower(varargin{6})
				case 'poisson'
					obj.distribution = 'Poisson';
					obj.a = [];
					obj.alpha = [];
					obj.R = [];
					obj.add = 0.0;
				case 'gaussian-independent'
					obj.distribution = 'Gaussian';
					obj.a = varargin{7}(1); % need checks
					obj.alpha = varargin{7}(2); % need checks
					obj.R = eye(obj.popSize);

					if length(varargin{7}) == 3
						obj.add = varargin{7}(3); % need checks
					else
						obj.add = 0.0;
					end
				case 'gaussian-uniform'
					obj.distribution = 'Gaussian';
					obj.a = varargin{7}(1); % need checks
					obj.alpha = varargin{7}(2); % need checks
					obj.R = varargin{7}(3) * ~eye(obj.popSize) + eye(obj.popSize);
					obj.add = 0.0;
				case 'gaussian-exponential'
					obj.distribution = 'Gaussian';
					obj.a = varargin{7}(1); % need checks
					obj.alpha = varargin{7}(2); % need checks
					c = varargin{7}(3); % need checks
					rho = varargin{7}(4); % need checks
					prefDiff = repmat(obj.preferredStimulus, 1, obj.popSize);
					prefDiff = prefDiff - prefDiff.';
					obj.R = c .* exp(-abs(double(prefDiff)) ./ rho) .* ~eye(obj.popSize) + eye(obj.popSize);
					obj.add = 0.0;
				case 'gaussian-gaussian'
					obj.distribution = 'Gaussian';
					obj.a = varargin{7}(1); % need checks
					obj.alpha = varargin{7}(2); % need checks
					c = varargin{7}(3); % need checks
					beta = 1.0 ./ degToRad(varargin{7}(4)).^2; % need checks
					prefDiff = repmat(obj.preferredStimulus, 1, obj.popSize);
					prefDiff = prefDiff - prefDiff.';
					obj.R = c .* exp((cosd(double(prefDiff)) - 1) .* beta) .* ~eye(obj.popSize) + eye(obj.popSize);
					obj.add = 0.0;
				case 'cercal'
					obj.distribution = 'Gaussian';
					obj.add = varargin{7}(1);
					obj.a = varargin{7}(2);
					obj.alpha = 0.5;
					obj.R = eye(obj.popSize);
					obj.exponent = 2.0;
				otherwise
					error([varargin{6} ' is not a valid variability regime'])
				end

			otherwise
				error('Wrong number of arguments')
			end
		end
		
		function varargout = mi(obj, method, stim, tol, maxiter)
			if sum(strcmp(method, {'quadrature' 'randMC' 'quasirandMC'})) == 0
				error([method ' is not a valid SSI calculation method'])
			end

			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(3) ' is not a SimulusEnsemble object'])
			end

			% obj.popSize x stim.n
			rMean = obj.integrationTime .* meanR(obj, stim);
			rMeanCell = squeeze(mat2cell(rMean, obj.popSize, ones(stim.n, 1)));

			% Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
			QCell1 = obj.Q(rMeanCell);

			% Compute Cholesky decomps of Q matrices
			cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);
			
			% Invert Q matrices
			invQCell = cellfun(@inv, QCell1, 'UniformOutput', false);
            clear QCell1
			cholInvQCell = cellfun(@chol, invQCell, 'UniformOutput', false);
            clear invQCell

			% Replicate cell arrays across stimulus ensemble
			%cholInvQCell = repmat(cholInvQCell', [1 stim.n]);
			%rMeanCella = repmat(rMeanCell', [1 stim.n]);
			
			% Define function for multivariate gaussian sampling
			% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
			
			if obj.truncate
                fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
            else
                fRand = @(m, c, z) m + c * z; % don't truncate
            end
			
			iter = 0; % iteration counter
            miEst = OnlineStats(1, maxiter);
            
            cpS = cumsum(stim.pS);
            
            cont = true;
			while cont
                iter = iter + 1;
                
                % Display progress every 100 iterations
				if ~mod(iter, 100)
					fprintf('mi()  iter: %d  val: %.4g  SEM: %.4g\n', iter, runMean, runSEM)
				end

				switch method
				case 'randMC'
					% Sample s from stimulus distribution
					[dummy, bin] = histc(rand(), cpS);
					bin = bin + 1;
					%s = double(stim.ensemble);
					%s = s(bin);

					% Sample r from response distribution
					% Generate vector of independent normal random numbers (mu=0, sigma=1)
					z = randn(obj.popSize, 1);
					% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
					% !!! NOTE NEGATIVE RESPONSES MAY BE TRUNCATED TO ZERO depending on value of obj.truncate !!!
					r = fRand(rMeanCell{bin}, cholQ{bin}, z);
				end

				% log P(r|s)
				% Replicate to form a stim.n x stim.n cell array of response vectors
				rCell = repmat({r}, [stim.n 1]);
				% Calculate response log probability densities
                lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCell, rMeanCell', cholInvQCell', {'inv'}));
                
				% log P(r,s')
				% Mutiply P(r|s) and P(s) to find joint distribution
				pS = stim.pS;
				lpRS = lpRgS + log(pS');

				% log P(r)
				% Calculate marginal by summing over s'
				lpR = logsumexp(lpRS);
				
				% log P(s)
				lpS = log(pS(bin));

				% MI in bits (convert from log_e to log_2)
				lpRS = lpRS(bin);
                
                % sample MI
                miEst.appendSample((lpRS - (lpR + lpS)) ./ log(2));
                
                % Test halting criteria (SEM, max iterations limit)
				cont = miEst.runSEM > tol & iter < maxiter;
                
                % Impose minimum iteration limit so we get a sensible estimate of SEM
                cont = cont | iter < 10;
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
            % Initialise wall clock if necessary
            try
                dummy = toc;
            catch
                tic
            end

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

			% Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
			QCell1 = obj.Q(rMeanCell);

			% Compute Cholesky decomps of Q matrices
			cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);

			% Invert Q matrices and compute Cholesky decomps
			invQCell = cellfun(@inv, QCell1, 'UniformOutput', false);
			cholInvQCell = cellfun(@chol, invQCell, 'UniformOutput', false);
			
			if ~isempty(n)
				% Create logical vector (mask) identifying neurons that are *not* part of the marginal SSI
				margMask = true(obj.popSize, 1);
				margMask(n) = false;
				
				% Get mean responses for each stimulus ordinate
				rMeanMargCell = cellfun(@(r) r(margMask), rMeanCell, 'UniformOutput', false);
				
				% Compute mean response dependent cov matrix stack Q
				QCellMarg1 = cellfun(@(q) q(margMask, margMask), QCell1, 'UniformOutput', false);
								
				% Invert Q matrices and compute Cholesky decomps
				invQCellMarg = cellfun(@inv, QCellMarg1, 'UniformOutput', false);
				cholInvQCellMarg = cellfun(@chol, invQCellMarg, 'UniformOutput', false);
			end
			
			% Define function for multivariate gaussian sampling
			% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
			% Comment as appropriate if you want to truncate at zero
			% This will mess up the Gaussianity
			
			if obj.truncate
                fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
            else
                fRand = @(m, c, z) m + c * z; % don't truncate
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
					fprintf('SSISS iter: %d of %d, SEM: %.4g\n', iter, maxiter, mean(Issi.runSEM))
				end

				switch method
				case 'randMC'
					% Sample r from response distribution
					% Generate vector of independent normal random numbers (mu=0, sigma=1)
					zCell = mat2cell(randn(obj.popSize, sMaskN), obj.popSize, ones(sMaskN, 1));
					% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
					% !!! NOTE NEGATIVE RESPONSES MAY BE TRUNCATED TO ZERO, SEE ABOVE !!!
					rCell = cellfun(fRand, rMeanCell(sMask), cholQ(sMask), zCell, 'UniformOutput', false); % stim.n cell array of obj.popSize vectors

				case 'quasirandMC'
					% 

				end

				% log P(r|s')
				% Calculate response probability densities
				lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCell, rMeanCell', cholInvQCell', {'inv'}));
                
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
					lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCellMarg, rMeanMargCell', cholInvQCellMarg', {'inv'}));
                    
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
                    runSEM = mean([Issi.runSEM IssiMarg.runSEM]);
                else
                    runSEM = mean(Issi.runSEM);
                end
                
                cont = all([runSEM > tol, iter < maxiter, toc < timeout]);
                
                % Impose minimum iteration limit so we get a valid estimate of SEM
                cont = cont | iter < 10;
			end

			fprintf('SSISS iter: %d  elapsed time: %.4f seconds\n', iter, toc)
            
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
			
			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(4) ' is not a SimulusEnsemble object'])
            end
            
            try
                dummy = toc;
            catch
                tic
            end
            
			switch method
			case 'analytic'	
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

				f = obj.integrationTime .* meanR(obj, stim);
				f_prime = obj.integrationTime .* dMeanR(obj, stim);

				g0 = f_prime ./ f;

				% ==============================================================
				% Correlation matrix, Covariance matrix, inverse and derivative.
				% Fourier tranform (eigenvalues)
				% ==============================================================

				fCell = squeeze(mat2cell(f, obj.popSize, ones(stim.n, 1)));
				f_primeCell = squeeze(mat2cell(f_prime, obj.popSize, ones(stim.n, 1)));
				g0Cell = squeeze(mat2cell(g0, obj.popSize, ones(stim.n, 1)));
				QCell1 = obj.Q(fCell);
				Q_inv = cellfun(@inv, QCell1, 'UniformOutput', false); % inverse
				k = repmat({obj.a(ones(obj.popSize,1))}, [1 stim.n]);
				k_prime = repmat({zeros(obj.popSize,1)}, [1, stim.n]);
				
				fQ_prime = @(q, kk, kp, gz) obj.alpha * (diag(gz) * q + q * diag(gz)) + (diag(kp ./ kk)) * q;
				Q_prime = cellfun(fQ_prime, QCell1, k, k_prime, g0Cell, 'UniformOutput', false); % derivative

				% ==============================================================
				% Fisher Information
				% I_mean, I_cov and approximation of I_cov
				% ==============================================================

				if obj.popSize > 1
					Info(1,:) = cellfun(@(fp, qi) fp' * qi * fp, f_primeCell, Q_inv);         % direct method
					%Info(2) = cellfun(@(fp, di) fp * di * fp', f_primeCell, D_inv);
				else
					Info(1,:) = cellfun(@(fp, q) fp.^2 / q, f_primeCell, QCell1);
				end

				Info(3,:) = cellfun(@(qp, qi) 0.5 * trace(qp * qi * qp * qi), Q_prime, Q_inv);  % direct method for Icov

				Info(5,:) = Info(1,:) + Info(3,:);                         % Fisher
                
				switch nargout
				case 1
					varargout = {Info(5,:)};
				case 2
					varargout = {Info(1,:) Info(3,:)};
				otherwise
					error('Wrong number of outputs')
				end

			case {'randMC' 'quasirandMC'}
                % Generate offset stimuli
                stim2 = stim;
                dTheta = 0.1 .* stim2.width;
                stim2.ensemble = stim2.ensemble + dTheta;
                
				% obj.popSize x stim.n
				rMean = obj.integrationTime .* obj.meanR(stim);
                rMeanCell = squeeze(mat2cell(rMean, obj.popSize, ones(stim.n, 1)));
                rMean2 = obj.integrationTime .* obj.meanR(stim2);
                rMeanCell2 = squeeze(mat2cell(rMean2, obj.popSize, ones(stim2.n, 1)));

				% Compute mean response dependent cov matrix stack Q [ (popSize x popSize) x stim.n ]
                QCell1 = obj.Q(rMeanCell);
                QCell2 = obj.Q(rMeanCell2);
                
                % Compute Cholesky decomps of Q matrices
                cholQ = cellfun(@(q) chol(q)', QCell1, 'UniformOutput', false);
                
                % Invert Q matrices and compute Cholesky decomps
                invQCell = cellfun(@inv, QCell1, 'UniformOutput', false);
                invQCell2 = cellfun(@inv, QCell2, 'UniformOutput', false);
                cholInvQCell = cellfun(@chol, invQCell, 'UniformOutput', false);
                cholInvQCell2 = cellfun(@chol, invQCell2, 'UniformOutput', false);
                
                % Concat (stim) and (stim + dTheta) means and covariances
                rMeanCellDiff = [rMeanCell ; rMeanCell2];
                cholInvQCellDiff = [cholInvQCell ; cholInvQCell2];
                
				% Define function for multivariate gaussian sampling
				% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
                
				if obj.truncate
				    fRand = @(m, c, z) max((m + c * z), 0.0); % truncate
				else
				    fRand = @(m, c, z) m + c * z; % don't truncate
                end
                
				iter = 0;
                cont = true;
                
                FI = OnlineStats(stim.n, maxiter);
                
                while cont
                    iter = iter + 1;
                    
                    if ~mod(iter, 100)
                        fprintf('Fisher iter: %d, SEM: %.4g\n', iter, mean(FI.runSEM))
                    end
                    
					switch method
					case 'randMC'
						% Sample r from response distribution
						% Generate vector of independent normal random numbers (mu=0, sigma=1)
						zCell = mat2cell(randn(obj.popSize, stim.n), obj.popSize, ones(stim.n, 1));
						% Multiply by Cholesky decomposition of cov matrix Q, and add in mean
						rCell = cellfun(fRand, rMeanCell, cholQ, zCell, 'UniformOutput', false); % stim.n cell array of {obj.popSize}
                        
					case 'quasirandMC'
						%
                    end

                    % log P(r|s)
                    % Calculate response probability densities
                    lpRgS = cell2mat(cellsxfun(@mvnormpdfln, rCell, rMeanCellDiff, cholInvQCellDiff, {'inv'})); % 1 x stim.n
                    
                    % d/ds log P(r|s)
					% Find derivative
                    dlpRgS = diff(lpRgS, 1, 1) ./ dTheta;

					% (d/ds log P(r|s))^2
                    FI.appendSample(dlpRgS .^ 2);
                    
					% Test halting criteria (SEM, max iterations limit, timeout)
                    runSEM = mean(FI.runSEM);
                    cont = all([runSEM > tol, iter < maxiter]);

                    % Impose minimum iteration limit so we get a sensible estimate of SEM
                    cont = cont | iter < 10;
                end
                
                fprintf('FI completed iter: %d\n', iter)
                
                % Trim sample array
                FI.trim;
                
                % Recalculate means cleanly
				varargout = {FI.mean};
			otherwise
				error([method ' is not a valid FI calculation method'])
			end
		end
		
		function varargout = margFisher(obj, nMarg, method, stim, tol, varargin)
			if isempty(varargin)
				option = 'raw';
			else
				option = varargin{1};
			end

			if ~isa(stim, 'StimulusEnsemble')
				error([inputname(4) ' is not a SimulusEnsemble object'])
			end

			% Compute FI of full population
			fullFisher = fisher(obj, method, stim, tol);
			
			% Create population of remaining neurons
			remainingNeurons = obj.remove(nMarg);
			
			% FI of remaining neurons
			remainderFisher = fisher(remainingNeurons, method, stim, tol);
			
			switch option
			case 'diff'
				varargout = {fullFisher - remainderFisher};
			case 'rootDiff'
				varargout = {sqrt(fullFisher - remainderFisher)};
			case 'diffRoot'
				varargout = {sqrt(fullFisher) - sqrt(remainderFisher)};
			case 'invDiffInv'
				varargout = {fullFisher .* remainderFisher ./ (remainderFisher - fullFisher)};
			case 'diffSD'
				varargout = {(fullFisher.^-0.5 - remainderFisher.^-0.5).^-2};
			case 'raw'
				varargout = {fullFisher remainderFisher};
			otherwise
				error('Invalid marginal fisher option')
			end
		end
		
		function ifish = Ifisher(obj, stim)
			fish = obj.fisher('analytic', stim, 0);
			pS = stim.pS;
			zOrds = find(fish == 0);
			
			if ~isempty(zOrds)
				fish(zOrds) = [];
				pS(zOrds) = [];
				pS = pS ./ sum(pS);
			end
			
            ifish = stim.entropy - 0.5 .* sum(pS .* (1.0 + log2(pi) + log2(exp(1)) - log2(fish)));
		end
		
		function ifish = mIfisher(obj, nMarg, stim)
			[fullFI remFI] = obj.margFisher(nMarg, 'analytic', stim, 0, 'raw');

            ifishFull = stim.entropy - 0.5 .* sum(stim.pS .* (1.0 + log2(pi) + log2(exp(1)) - log2(fullFI)));
            ifishRem = stim.entropy - 0.5 .* sum(stim.pS .* (1.0 + log2(pi) + log2(exp(1)) - log2(remFI)));

			ifish = ifishFull - ifishRem;
		end
		
		function ssif = SSIfisher(obj, n, fisherMethod, stim, tol)
		
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
				[fullFI remFI] = obj.margFisher(n, fisherMethod, stim, tol, 'raw');
				
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
		
		function obj = gainadapt(obj, width, amnt, centre)
			obj.maxRate = obj.maxRate .* (1 - amnt .* exp(-(1.0 ./ degToRad(width).^2) .* (1 - cosd(double(obj.preferredStimulus - centre)))));
		end
		
		function p = pOfR(obj, response, stimulusSet)
			if ~(isnumeric(response) && isvector(response) && length(response) == obj.popSize)
				error([inputname(2) ' is not a valid response'])
			end

			if ~(isnumeric(stimulusSet) && size(stimulusSet, 1) == obj.dimensionality)
				error([inputname(3) ' is not a valid stimulus set'])
			end

			r = repmat(reshape(response, length(response), 1), 1, length(stimulusSet));
			rMean = meanR(obj, stimulusSet);
			r = r - rMean;
			
			prob = zeros(1,length(stimulusSet));
			
			for s = 1 : length(stimulusSet)
				Q = obj.varA * obj.varR .* (rMean(:,s) * rMean(:,s)') .^ obj.varAlpha;
				normFactor = 1.0 / ((2.0 * pi)^(obj.popSize / 2.0) * det(Q)^0.5);
				prob(1,s) = normFactor * exp(r(:,s)' * (Q \ r(:,s)));
			end

			p = prob;
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
			maxRateStr = char(obj.maxRate);
			backgroundRateStr = char(obj.backgroundRate);
			rStr = display(obj.R);

			formatStr = ['Population size: %.0f\n'...
						 'Stimulus dimensionality: %.0f\n'...
						 'Preferred stimuli:\n'...
						  prefStr '\n'...
						 'Maximum firing rate (Hz):\n'...
						  maxRateStr '\n'...
						 'Background firing rate (Hz):\n'...
						  backgroundRateStr '\n'...
						 'Integration time: %.3f s\n'...
						 'Variability distribution: %s\n'...
						 'Variability coefficient a: %.3f\n'...
						 'Variability exponent alpha: %.3f\n'...
						 'Correlation matrix:\n'];

			retStr = strvcat(sprintf(formatStr, obj.popSize, double(obj.dimensionality), obj.integrationTime, obj.distribution, obj.a, obj.alpha), rStr);
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
			
			if length(obj.maxRate) > 1
				obj.maxRate = obj.maxRate(margMask);
			end
			
			if length(obj.backgroundRate) > 1
				obj.backgroundRate = obj.backgroundRate(margMask);
			end
			
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

end
