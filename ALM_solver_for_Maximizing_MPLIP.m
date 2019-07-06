function [x, fval] = ALM_solver_for_Maximizing_MPLIP(amp, phase, ...
    relFreqConst,numSources,xcoords, ycoords, zcoords, sigma, ...
    epsilon, alpha, beta, maxIter, lb, ub)
% Maximizing Monochromatic Polarized Light Interference Patterns Using
% Augmented Lagrangian Method.

%#codegen

%% Optimization problem
% Construct the objective function
origins = [xcoords ycoords zcoords];

    function f = objFcn(x)
        f = - waveIntensity(x, amp,phase, ...
            relFreqConst,numSources,origins);
    end

% Construct the inequality constraints
ineqConsNum = 7;
% ineqCons = cell(ineqConsNum, 1);

g1 = @(x,y)  (x-xcoords(1)).^2 + (y-ycoords(1)).^2 - 9;
g2 = @(x,y)  (x-xcoords(2)).^2 + (y-ycoords(2)).^2 - 9;
g3 = @(x,y)  (x-xcoords(3)).^2 + (y-ycoords(3)).^2 - 9;
% ineqCons{1} = @(x)-g1(x(1), x(2));
% ineqCons{2} = @(x)-g2(x(1), x(2));
% ineqCons{3} = @(x)-g3(x(1), x(2));
%
% ineqCons{4} = @(x)(x(1) - lb(1));
% ineqCons{5} = @(x)(ub(1) - x(1));
% ineqCons{6} = @(x)(x(2) - lb(2));
% ineqCons{7} = @(x)(ub(2) - x(2));
ineqCons = {@(x)-g1(x(1), x(2)),...
    @(x)-g2(x(1), x(2)), ...
    @(x)-g3(x(1), x(2)), ...
    @(x)(x(1) - lb(1)), ...
    @(x)(ub(1) - x(1)), ...
    @(x)(x(2) - lb(2)), ...
    @(x)(ub(2) - x(2))};

% Construct the equality constraints
eqConsNum = 0;
eqCons = {};

nvars = 2;
w = ones(ineqConsNum, 1).*1.3;
v = ones(eqConsNum, 1).*1.3;

%% ALM algorithm kernel
% Step 1
% Initialize the epoch flag to 1.
x0 = [0.5; 0.5];
x = x0;
fval =  0;
k = 1;
x_k_1 = x0;
w_k = w;
v_k = v;
sigma_k = sigma;

% Step 2
while (k <=maxIter)
    % Construct the Augmented Lagrangian Function as the objective function of
    % a unconstrained optimization problem.
    [x_k, ~, ~] = PSO_alg_solver(w_k, v_k, sigma_k, nvars, lb, ub, 'off');
    
    % Step 3
    % Check the stop condition
    norm_h_k = norm_eqCons(x_k, eqCons, eqConsNum);
    if norm_h_k < epsilon
        x = [x_k(1);x_k(2)];
        fval = objFcn(x_k);
        break;
    end
    
    % Step 4
    % Adapt the sigma
    % Check the stop condition
    norm_h_k_1 = norm_eqCons(x_k_1, eqCons, eqConsNum);
    if norm_h_k >= beta*norm_h_k_1
        sigma_k = alpha.*sigma_k;
    end
    
    % Step 5
    % Update the multipliers
    [w_new, v_new] = multipliers_update(x_k, w_k, v_k, ...
        eqCons, eqConsNum, ineqCons, ineqConsNum, sigma_k);
    w_k = w_new;
    v_k = v_new;
    
    x_k_1 = [x_k(1);x_k(2)];
    
    k = k + 1;
end

%% Formulation of augmented Lagrangian function
    function f = Lagrange_Fcn(x, w, v, sigma)
        G = zeros(7, 1, 'double');
        H = zeros(0, 1, 'double');
        for i = 1:7
            G(i) = (max([0, w(i) - sigma.*ineqCons{i}(x)])).^2 - w(i).^2;
        end
        for j = 1:0
            H(j) = eqCons{j}(x);
        end
        f = objFcn(x) + (1./(2*sigma)).*(sum(G)) ...
            - sum(v.*H) + (sigma./2).*sum(H.^2);
    end

%% PSO solver for unconstrained problem
    function [x, fval, exitFlag] = PSO_alg_solver(w_k, v_k, sigma_k, ...
            nvars, lb, ub, Display)
        % Construct augmented Lagrangian function.
        objFcn = @(x)Lagrange_Fcn(x, w_k, v_k, sigma_k);
        
%         output = struct('rngstate', rng, 'iterations', 0, 'funccount', 0, 'message', ' ');
%         output = struct('rngstate', rng, 'iterations', 0, 'funccount', 0);
        [x,fval] = deal([]);
        
        options = struct;
        % options.Display = 'iter';
        options.DisplayInterval = 1;
        % switch  options.Display
        switch  Display
            case {'off','none'}
                options.Verbosity = 0;
            case 'final'
                options.Verbosity = 1;
            case 'iter'
                options.Verbosity = 2;
            otherwise
                options.Verbosity = 0;
        end
        
        options.InertiaRange = [0.1000 1.1000];
        options.InitialSwarm = [];
        options.SwarmSize = min(500, 100*nvars);
        options.InitialSwarmSpan = 2000;
        options.MaxIter = 200*nvars;
        options.MaxTime = Inf;
        options.MinFractionNeighbors = 0.2500;
        options.ObjectiveLimit = -Inf;
        options.SelfAdjustment = 1.4900;
        options.SocialAdjustment = 1.4900;
        options.StallIterLimit = 20;
        options.StallTimeLimit = Inf;
        options.TolFun = 1.0000e-6;
        options.TolFunValue = 1.0000e-6;
        
        % Args lb & ub
        [lbColumn, ubColumn, msg, exitFlag] = checkbound(lb, ub, nvars);
        % checkbound returns column vectors, but we to convert that to a row vector
        % and repeat it across SwarmSize rows.
        lbRow = lbColumn';
        ubRow = ubColumn';
        if exitFlag < 0
            % output.message = msg; % This is the first message, so no need to append
            if options.Verbosity > 0
                fprintf('%s\n', msg);
            end
            return
        end
        
        % Perform check on initial population, fvals, and range
        options = initialParticleCheck(options);
        
        exitFlag = 0;
        
        % Get algorithmic options
        numParticles = options.SwarmSize;
        cSelf = options.SelfAdjustment;
        cSocial = options.SocialAdjustment;
        minNeighborhoodSize = max(2,floor(numParticles*options.MinFractionNeighbors));
        minInertia = options.InertiaRange(1);
        maxInertia = options.InertiaRange(2);
        lbMatrix = repmat(lbRow, numParticles, 1);
        ubMatrix = repmat(ubRow, numParticles, 1);
        
        
        % Create initial state: particle positions & velocities, fvals, status data
        state = makeState(nvars,lbMatrix,ubMatrix,objFcn,options);
        bestFvals = min(state.Fvals);
        % Create a vector to store the last StallIterLimit bestFvals.
        % bestFvalsWindow is a circular buffer, so that the value from the i'th
        % iteration is stored in element with index mod(i-1,StallIterLimit)+1.
        bestFvalsWindow = nan(options.StallIterLimit, 1);
        
        % Initialize adaptive parameters:
        %   initial inertia = maximum *magnitude* inertia
        %   initial neighborhood size = minimum neighborhood size
        adaptiveInertiaCounter = 0;
        if all(options.InertiaRange >= 0)
            adaptiveInertia = maxInertia;
        elseif all(options.InertiaRange <= 0)
            adaptiveInertia = minInertia;
        else
            % checkfield should prevent InertiaRange from having positive and
            % negative vlaues.
            assert(false, 'PSO:particleswarm:invalidInertiaRange', ...
                'The InertiaRange option should not contain both positive and negative numbers.');
        end
        adaptiveNeighborhoodSize = minNeighborhoodSize;
        
        % Setup display header
        if  options.Verbosity > 1
            fprintf('\n                                 Best            Mean     Stall\n');
            fprintf(  'Iteration     f-count            f(x)            f(x)    Iterations\n');
            fprintf('%5.0f         %7.0f    %12.4g    %12.4g    %5.0f\n', ...
                0, state.FunEval, bestFvals, mean(state.Fvals), 0);
        end
        
        % Run the main loop until some exit condition becomes true
%         message = '';
        while (~exitFlag)
            state.Iteration = state.Iteration + 1;
            
            % Generate a random neighborhood for each particle that includes
            % the particle itself
            neighborIndex = zeros(numParticles, adaptiveNeighborhoodSize);
            neighborIndex(:, 1) = 1:numParticles; % First neighbor is self
            for i = 1:numParticles
                % Determine random neighbors that exclude the particle itself,
                % which is (numParticles-1) particles
                neighbors = randperm(numParticles-1, adaptiveNeighborhoodSize-1);
                % Add 1 to indicies that are >= current particle index
                iShift = neighbors >= i;
                neighbors(iShift) = neighbors(iShift) + 1;
                neighborIndex(i,2:end) = neighbors;
            end
            % Identify the best neighbor
            [~, bestRowIndex] = min(state.IndividualBestFvals(neighborIndex), [], 2);
            % Create the linear index into neighborIndex
            bestLinearIndex = (bestRowIndex.'-1).*numParticles + (1:numParticles);
            bestNeighborIndex = neighborIndex(bestLinearIndex);
            randSelf = rand(numParticles, nvars);
            randSocial = rand(numParticles, nvars);
            
            % Note that velocities and positions can become infinite if the
            % inertia range is too large or if the objective function is badly
            % behaved.
            
            % Update the velocities
            newVelocities = adaptiveInertia*state.Velocities + ...
                cSelf*randSelf.*(state.IndividualBestPositions-state.Positions) + ...
                cSocial*randSocial.*(state.IndividualBestPositions(bestNeighborIndex, :)-state.Positions);
            tfValid = all(isfinite(newVelocities), 2);
            state.Velocities(tfValid,:) = newVelocities(tfValid,:);
            % Update the positions
            newPopulation = state.Positions + state.Velocities;
            tfInvalid = ~isfinite(newPopulation);
            newPopulation(tfInvalid) = state.Positions(tfInvalid);
            % Enforce bounds, setting the corresponding velocity component to
            % zero if a particle encounters a lower/upper bound
            tfInvalid = newPopulation < lbMatrix;
            if any(tfInvalid(:))
                newPopulation(tfInvalid) = lbMatrix(tfInvalid);
                state.Velocities(tfInvalid) = 0;
            end
            tfInvalid = newPopulation > ubMatrix;
            if any(tfInvalid(:))
                newPopulation(tfInvalid) = ubMatrix(tfInvalid);
                state.Velocities(tfInvalid) = 0;
            end
            state.Positions = newPopulation;
            
            % Objective evaluation
            state.Fvals = fcnvectorizer(objFcn, state.Positions);
            state.FunEval = state.FunEval + numParticles;
            
            % Remember the best fvals and positions
            tfImproved = state.Fvals < state.IndividualBestFvals;
            state.IndividualBestFvals(tfImproved) = state.Fvals(tfImproved);
            state.IndividualBestPositions(tfImproved, :) = state.Positions(tfImproved, :);
            bestFvalsWindow(1+mod(state.Iteration-1,options.StallIterLimit)) = min(state.IndividualBestFvals);
            
            % Keep track of improvement in bestFvals and update the adaptive
            % parameters according to the approach described in S. Iadevaia et
            % al. Cancer Res 2010;70:6704-6714 and M. Liu, D. Shin, and H. I.
            % Kang. International Conference on Information, Communications and
            % Signal Processing 2009:1-5.
            
            newBest = min(state.IndividualBestFvals);
            if isfinite(newBest) && newBest < bestFvals
                bestFvals = newBest;
                state.LastImprovement = state.Iteration;
                %state.LastImprovementTime = toc(state.StartTime);
                adaptiveInertiaCounter = max(0, adaptiveInertiaCounter-1);
                adaptiveNeighborhoodSize = minNeighborhoodSize;
            else
                adaptiveInertiaCounter = adaptiveInertiaCounter+1;
                adaptiveNeighborhoodSize = min(numParticles, adaptiveNeighborhoodSize+minNeighborhoodSize);
            end
            
            % Update the inertia coefficient, enforcing limits (Since inertia
            % can be negative, enforcing both upper *and* lower bounds after
            % multiplying.)
            if adaptiveInertiaCounter < 2
                adaptiveInertia = max(minInertia, min(maxInertia, 2*adaptiveInertia));
            elseif adaptiveInertiaCounter > 5
                adaptiveInertia = max(minInertia, min(maxInertia, 0.5*adaptiveInertia));
            end
            
            % check to see if any stopping criteria have been met
%             [exitFlag, message] = stopParticleswarm(options,state,bestFvalsWindow);
            [exitFlag, ~] = stopParticleswarm(options,state,bestFvalsWindow);
%             output.message = message;
        end % End while loop
        
        % Find and return the best solution
        [fval,indexBestFval] = min(state.IndividualBestFvals);
        x = state.IndividualBestPositions(indexBestFval,:);
        
        % Update output structure
%         iterations = state.Iteration;
%         funccount   = state.FunEval;
%         output = struct('rngstate', rng, ...
%             'iterations', iterations, ...
%             'funccount', funccount, ...
%             'message', message);
    end

    function [lb,ub,msg,exitflag] = checkbound(lbin,ubin,nvars)
        %CHECKBOUND Move the initial point within the (valid) bounds.
        %   [X,LB,UB,X,FLAG] = CHECKBOUNDS(X0,LB,UB,nvars) checks that the upper and lower
        %   bounds are valid (LB <= UB) and the same length as X (pad with -inf/inf
        %   if necessary); warn if too long.  Also make LB and UB vectors if not
        %   already. Finally, inf in LB or -inf in UB throws an error.
        
        msg = 'Success';
        exitflag = 1;
        
        % Turn into column vectors
        lb = lbin(:);
        ub = ubin(:);
        lenlb = length(lb);
        lenub = length(ub);
        
        % Check maximum length
        if lenlb > nvars
            warning(message('PSO:checkbound:lengthOfLowerBound'));
            lb = lb(1:nvars);
            lenlb = nvars;
        elseif lenlb < nvars
            lb = [lb; -inf*ones(nvars-lenlb,1)];
            lenlb = nvars;
        end
        
        if lenub > nvars
            warning(message('PSO:checkbound:lengthOfUpperBound'));
            ub = ub(1:nvars);
            lenub = nvars;
        elseif lenub < nvars
            ub = [ub; inf*ones(nvars-lenub,1)];
            lenub = nvars;
        end
        
        % Check feasibility of bounds
        len = min(lenlb,lenub);
        if any( lb( (1:len)' ) > ub( (1:len)' ) )
            count = full(sum(lb>ub));
            count = int64(count);
            if count == 1
                msg=sprintf(['Exiting due to infeasibility:  %i lower bound exceeds the' ...
                    ' corresponding upper bound.'],count);
            else
                msg=sprintf(['Exiting due to infeasibility:  %i lower bounds exceed the' ...
                    ' corresponding upper bounds.'],count);
            end
            exitflag = -2;
        end
        % check if -inf in ub or inf in lb
        if any(eq(ub, -inf))
            error(message('PSO:checkbound:infUpperBound'));
        elseif any(eq(lb,inf))
            error(message('PSO:checkbound:infLowerBound'));
        end
    end % End of checkbound

    function options = initialParticleCheck(options)
        % checkfield already checks that InitialSwarm is a matrix with nvars
        % columns
        
        numInitPositions = size(options.InitialSwarm, 1);
        
        % No tests if initial positions is empty
        if numInitPositions == 0
            return
        end
        
        % Warn if too many positions were specified.
        numParticles  = options.SwarmSize;
        if numInitPositions > numParticles
            warning(message('PSO:particleswarm:initialSwarmLength'));
            options.InitialSwarm(numParticles+1:numInitPositions,:) = [];
        end
        
    end % End of initialParticleCheck


    function state = makeState(nvars, lbMatrix, ubMatrix, objFcn, options)
        % Create an initial set of particles and objective function values
        
        % makeState needs the vector of bounds, not the expanded matrix.
        lb_local = lbMatrix(1,:);
        ub_local = ubMatrix(1,:);
        
        % A variety of data used in various places
        state = struct;
        state.Iteration = 0; % current generation counter
        %state.StartTime = tic; % tic identifier
        state.StopFlag = false; % OutputFcns flag to end the optimization
        state.LastImprovement = 1; % generation stall counter
        state.LastImprovementTime = 0; % stall time counter
        state.FunEval = 0;
        numParticles = options.SwarmSize;
        
        % If InitialSwarm is partly empty use the creation function to generate
        % population (CreationFcn can utilize InitialSwarm)
        if numParticles ~= size(options.InitialSwarm,1)
            problemStruct = struct;
            problemStruct.rngstate = [];
            problemStruct.solver = 'particleswarm';
            problemStruct.objective = objFcn;
            problemStruct.lb = lb_local;
            problemStruct.ub = ub_local;
            problemStruct.nvars = nvars;
            problemStruct.options = options;
            
            Positions = pswcreationuniform(problemStruct);
        else % the initial swarm was passed in!
            Positions = options.InitialSwarm;
        end
        
        % Enforce bounds
        if any(any(Positions < lbMatrix)) || any(any(Positions > ubMatrix))
            Positions = max(lbMatrix, Positions);
            Positions = min(ubMatrix, Positions);
            if options.Verbosity > 1
                fprintf(getString(message('PSO:particleswarm:shiftX0ToBnds')));
            end
        end
        
        % Initialize velocities by randomly sampling over the smaller of
        % options.InitialSwarmSpan or ub-lb. Note that min will be
        % InitialSwarmSpan if either lb or ub is not finite.
        vmax = min(ub_local-lb_local, options.InitialSwarmSpan);
        Velocities = repmat(-vmax,numParticles,1) + ...
            repmat(2*vmax,numParticles,1) .* rand(numParticles,nvars);
        
        % Calculate the objective function for all particles.
        firstFval = objFcn(Positions(1,:));
        % User-provided objective function should return a scalar
        if numel(firstFval) ~= 1
            error(message('PSO:particleswarm:objectiveCheck'));
        end
        fvals = fcnvectorizer(objFcn, Positions(2:end,:));
        % Concatenate the fvals of the first particle to the rest
        Fvals = [firstFval; fvals];
        
        state.FunEval = numParticles;
        
        state.IndividualBestFvals = Fvals;
        state.IndividualBestPositions = Positions;
        state.Velocities = Velocities;
        state.Positions = Positions;
        state.Fvals = Fvals;
    end % End of makeState

    function swarm = pswcreationuniform(problemStruct)
        %PSWCREATIONUNIFORM Creates the initial positions for PARTICLESWARM algorithm.
        
        nvars_local = problemStruct.nvars;
        options = problemStruct.options;
        % Determine finite bounds for the initial particles based on the problem's
        % bounds and options.InitialSwarmSpan.
        [lb_local,ub_local] = determinePositionInitBounds(problemStruct.lb, problemStruct.ub, ...
            options.InitialSwarmSpan);
        
        numParticles = options.SwarmSize;
        numInitPositions = size(options.InitialSwarm, 1);
        numPositionsToCreate = numParticles - numInitPositions;
        
        % Initialize particles to be created
        swarm = zeros(numParticles,nvars_local);
        
        % Use initial particles provided already
        if numInitPositions > 0
            swarm(1:numInitPositions,:) = options.InitialSwarm;
        end
        
        % Create remaining particles, randomly sampling within lb and ub
        span = ub_local - lb_local;
        swarm(numInitPositions+1:end,:) = repmat(lb_local,numPositionsToCreate,1) + ...
            repmat(span,numPositionsToCreate,1) .* rand(numPositionsToCreate,nvars_local);
        
        % Error if any values are not finite
        if ~all(isfinite(swarm(:)))
            error(message('PSO:pswcreationuniform:positionNotFinite'));
        end
    end

    function [lb,ub] = determinePositionInitBounds(lb,ub,initialSwarmSpan)
        % Update lb and ub using positionInitSpan, so that initial bounds are
        % always finite
        lbFinite = isfinite(lb);
        ubFinite = isfinite(ub);
        lbInf = ~lbFinite;
        ubInf = ~ubFinite;
        
        % If lb and ub are both finite, do not update the bounds.
        
        % If lb & ub are both infinite, center the range around 0.
        idx = lbInf & ubInf;
        lb(idx) = -initialSwarmSpan(idx)/2;
        ub(idx) = initialSwarmSpan(idx)/2;
        
        % If only lb is finite, start the range at lb.
        idx = lbFinite & ubInf;
        ub(idx) = lb(idx) + initialSwarmSpan(idx);
        
        % If only ub is finite, end the range at ub.
        idx = lbInf & ubFinite;
        lb(idx) = ub(idx) - initialSwarmSpan(idx);
    end

    function fvals = fcnvectorizer(objfcn, X)
        particle_num = size(X, 1);
        fvals = zeros(particle_num,1);
        for i = 1:particle_num
            fvals(i) = objfcn(X(i,:));
        end
    end

    function [exitFlag,reasonToStop] = stopParticleswarm(options,state,bestFvalsWindow)
        iteration = state.Iteration;
        
        iterationIndex = 1+mod(iteration-1,options.StallIterLimit);
        bestFval = bestFvalsWindow(iterationIndex);
        if options.Verbosity > 1 && ...
                mod(iteration,options.DisplayInterval)==0 && ...
                iteration > 0
            FunEval  = state.FunEval;
            MeanFval = meanf(state.Fvals);
            StallGen = iteration  - state.LastImprovement;
            fprintf('%5.0f         %7.0f    %12.4g    %12.4g    %5.0f\n', ...
                iteration, FunEval, bestFval, MeanFval, StallGen);
        end
        
        % Compute change in fval and individuals in last 'Window' iterations
        Window = options.StallIterLimit;
        if iteration > Window
            % The smallest fval in the window should be bestFval.
            % The largest fval in the window should be the oldest one in the
            % window. This value is at iterationIndex+1 (or 1).
            if iterationIndex == Window
                % The window runs from index 1:iterationIndex
                maxBestFvalsWindow = bestFvalsWindow(1);
            else
                % The window runs from [iterationIndex+1:end, 1:iterationIndex]
                maxBestFvalsWindow = bestFvalsWindow(iterationIndex+1);
            end
            funChange = abs(maxBestFvalsWindow-bestFval)/max(1,abs(bestFval));
        else
            funChange = Inf;
        end
        
        reasonToStop = '';
        exitFlag = 0;
        if state.Iteration >= options.MaxIter
            reasonToStop = getString(message('PSO:particleswarm:ExitMaxIter'));
            exitFlag = 2;
%         elseif toc(state.StartTime) > options.MaxTime
%             reasonToStop = getString(message('PSO:particleswarm:ExitMaxTime'));
%             exitFlag = -5;
%         elseif (toc(state.StartTime)-state.LastImprovementTime) > options.StallTimeLimit
%             reasonToStop = getString(message('PSO:particleswarm:ExitStallTimeLimit'));
%             exitFlag = -4;
        elseif bestFval <= options.ObjectiveLimit
            reasonToStop = getString(message('PSO:particleswarm:ExitObjectiveLimit'));
            exitFlag = -3;
        elseif state.StopFlag
            reasonToStop = getString(message('PSO:particleswarm:ExitOutputPlotFcn'));
            exitFlag = -1;
        elseif funChange <= options.TolFunValue
            reasonToStop = getString(message('PSO:particleswarm:ExitTolFun'));
            exitFlag = 1;
        end
        
        if ~isempty(reasonToStop) && options.Verbosity > 0
            fprintf('%s\n',reasonToStop);
            return
        end
        
        % Print header again
        if options.Verbosity > 1 && rem(iteration,30*options.DisplayInterval)==0 && iteration > 0
            fprintf('\n                                 Best            Mean     Stall\n');
            fprintf('Iteration     f-count            f(x)            f(x)    Iterations\n');
        end
        
    end

    function m = meanf(x)
        tfValid = ~isnan(x);
        n = sum(tfValid);
        if n==0
            % prevent divideByZero warnings
            m = NaN;
        else
            % Sum up non-NaNs, and divide by the number of non-NaNs.
            m = sum(x(tfValid)) ./ n;
        end
    end % End of StopParticleswarm

    function msg_str = getString(msg_obj)
        %msg_str = msg_obj.Identifier;
        msg_str = msg_obj;
    end

end

function [w_new, v_new] = multipliers_update(x_k, w_k, v_k, ...
    eqCons, eqConsNum, ineqCons, ineqConsNum, sigma)
% Allocate memory
w_new = w_k;
v_new = v_k;
for i = 1:7
    w_new(i) = max([0, w_k(i) - sigma.*ineqCons{i}(x_k)]);
end
for j = 1:0
    v_new(j) = v_k(j) - sigma.*eqCons{j}(x_k);
end
end

function norm_h = norm_eqCons(x_k, eqCons, eqConsNum)
H = zeros(0, 1, 'double');
for i = 1:0
    H(i) = eqCons{i}(x_k);
end
norm_h = norm(H);
end

function str_o = message(str)
str_o = str;
end