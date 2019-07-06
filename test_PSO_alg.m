%% Test PSO Algorithm Prototype

clear 
close all
clc

%% A simple testing optimization problem
% The following can be compared with the built-in particleswarm solver
% provided by MATLAB.
%      objFcn = @dejong5fcn;
%      nvars = 2;
%      lb = [-64 -64];
%      ub = [64 64];
%      [x,fval] = particleswarm(objFcn,nvars,lb,ub)
objFcn = @dejong5fcn;
nvars = 2;
lb = [-64 -64];
ub = [64 64];
[x, fval, exitFlag, output] = PSO_alg(objFcn, nvars, lb, ub, 'iter')

%% A more complex testing optimization problem
%
% This example shows how to locate the maximum
% of an electromagnetic interference pattern. For simplicity of modeling,
% the pattern arises from monochromatic polarized light spreading out from
% point sources.
% 
% The electric field due to source i measured in the direction of
% polarization at point x and time t is 
%
% $$E_i = \frac{A_i}{d_i(x)} \sin(\phi_i + \omega (t - d_i(x)/c) ),$$ 
%
% where $\phi_i$ is the phase at time zero for source $i$, $c$ is the speed
% of light, $\omega$ is the frequency of the light, $A_i$ is the amplitude
% of source $i$, and $d_i(x)$ is the distance from source $i$ to $x$.
%  
% For a fixed point $x$ the intensity of the light is the time average of
% the square of the net electric field. The net electric field is sum of
% the electric fields due to all sources. The time average depends only on
% the sizes and relative phases of the electric fields at $x$. To calculate
% the net electric field, add up the individual contributions using the
% phasor method. For phasors, each source contributes a vector. The length
% of the vector is the amplitude divided by distance from the source, and
% the angle of the vector, $\phi_i - \omega d_i(x)/c$ is the phase at the
% point.
%
% For this example, we define three point sources with the same frequency
% ($\omega$) and amplitude ($A$), but varied initial phase ($\phi_i$). We
% arrange these sources on a fixed plane.

% Frequency is proportional to the number of peaks
relFreqConst = 2*pi*2.5;
amp = 2.2;
phase = -[0; 0.54; 2.07];

numSources = 3;
height = 3;

% All point sources are aligned at [x_i,y_i,z]
xcoords = [2.4112
           0.2064
           1.6787];
ycoords = [0.3957
           0.3927
           0.9877];
zcoords = height*ones(numSources,1);          

origins = [xcoords ycoords zcoords];

%% Visualize the Interference Pattern
% Now let's visualize a slice of the interference pattern on the plane z = 0.
%
% As you can see from the plot below, there are many peaks and valleys
% indicating constructive and destructive interference.

% Pass additional parameters via an anonymous function:
waveIntensity_x = @(x) waveIntensity(x,amp,phase, ...
    relFreqConst,numSources,origins);
% Generate the grid
[X,Y] = meshgrid(-4:0.035:4,-4:0.035:4);
% Compute the intensity over the grid
Z = arrayfun(@(x,y) waveIntensity_x([x y]),X,Y);
% Plot the surface and the contours
figure
surf(X,Y,Z,'EdgeColor','none')
xlabel('x')
ylabel('y')
zlabel('intensity')

%% Posing the Optimization Problem
% We are interested in the location where this wave intensity reaches
% its highest peak. 
%
% The wave intensity ($I$) falls off as we move away from the source 
% proportional to $1/d_i(x)$. Therefore, let's restrict the space of
% viable solutions by adding constraints to the problem.
%
% If we limit the exposure of the sources with an aperture, then we can
% expect the maximum to lie in the intersection of the projection of the
% apertures onto our observation plane. We model the effect of an aperture
% by restricting the search to a circular region centered at each source. 
%
% We also restrict the solution space by adding bounds to the problem.
% Although these bounds may be redundant (given the nonlinear constraints),
% they are useful since they restrict the range in which start points are
% generated.
%
% Now our problem has become:
%
% $$ \max_{x,y} I(x,y) $$
%
% subject to
%
% $$ (x - x_{c1})^2 + (y - y_{c1})^2 \le r_1^2 $$
%
% $$ (x - x_{c2})^2 + (y - y_{c2})^2 \le r_2^2 $$
%
% $$ (x - x_{c3})^2 + (y - y_{c3})^2 \le r_3^2 $$
%
% $$-0.5 \leq x \leq 3.5$$
%
% $$-2 \leq y \leq 3$$
%
% where $(x_{cn},y_{cn})$ and $r_n$ are the coordinates and aperture radius
% of the $n^{th}$ point source, respectively. Each source is given an
% aperture with radius 3. The given bounds encompass the feasible region.
%
% The objective ($I(x,y)$) and nonlinear constraint functions are defined
% in separate MATLAB(R) files, |waveIntensity.m| and |apertureConstraint.m|,
% respectively, which are listed at the end of this example.
%

%% Visualization the objective function
% Now let's visualize the contours of our interference pattern with the
% nonlinear constraint boundaries superimposed. The feasible region is the
% interior of the intersection of the three circles (yellow, green, and
% blue). The bounds on the variables are indicated by the dashed-line box.
%

% Visualize the contours of our interference surface
domain = [-3 5.5 -4 5];
figure;
ezcontour(@(X,Y) arrayfun(@(x,y) waveIntensity_x([x y]),X,Y),domain,150);
hold on
title('Pattern Contours of The Interference Surface')

%% Setting Up and Solving the Problem with a Local Solver
% Given the nonlinear constraints, we need a constrained nonlinear solver,
% namely, |fmincon|. 
%
% Let's set up a problem structure describing our optimization problem. We
% want to maximize the intensity function, so we negate the values returned
% form |waveIntensity|. Let's choose an arbitrary start point that happens
% to be near the feasible region. 
%
% For this small problem, we'll use |fmincon|'s SQP algorithm.

% Pass additional parameters via an anonymous function:
apertureConstraint_x = @(x) apertureConstraint(x,xcoords,ycoords);

% Set up fmincon's options
x0 = [3 -1];
opts = optimoptions('fmincon','Algorithm','sqp');
problem = createOptimProblem('fmincon','objective', ...
    @(x) -waveIntensity_x(x),'x0',x0,'lb',lb,'ub',ub, ...
    'nonlcon',apertureConstraint_x,'options',opts);

% Call fmincon
[xlocal,fvallocal] = fmincon(problem)
%% 
% Now, let's see how we did by showing the result of |fmincon| in our contour
% plot. Notice that |fmincon| did not reach the global maximum, which is also
% annotated on the plot. Note that we'll only plot the bound that was
% active at the solution.

[~,maxIdx] = max(Z(:));
xmax = [X(maxIdx),Y(maxIdx)];
figure
contour(X,Y,Z)
hold on

% Show bounds
line([lb(1) lb(1)],[lb(2) ub(2)],'LineStyle','--')

% Plot PSO results using a circle marker
P1 = scatter(xlocal(1), xlocal(2), 's', ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor',[0 0 1],'LineWidth',1.1);

% Create textarrow showing the location of xglobal
annotation('textarrow',[0.44 0.50],[0.63 0.58],'TextEdgeColor',[0 0 0],...
    'TextBackgroundColor',[1 1 1],'FontSize',12,'String',{'Global Max'});
legend(P1, 'MATLAB:fmincon()', 'Location', 'best');
axis([-1 3.75 -3 3])
%% Using |PSO algorithm|
% Given an arbitrary initial guess, |fmincon| gets stuck at a nearby local
% maximum. Finding tight bounds can be difficult to do in practice, when not much is
% known about the objective function or constraints. In general though, we
% may be able to guess a reasonable region in which we would like to
% restrict the set of start points.

rng(4,'twister') % for reproducibility

% Run GlobalSearch
nvars = 2;
lb = -5*ones(2,1);
ub = 5*ones(2,1);
objFcn = @(x) -waveIntensity_x(x);
tic;
[x, fval, exitFlag, output] = PSO_alg(objFcn, nvars, lb, ub, 'iter')
toc


%% Examining Results
% Show the contours
figure
contour(X,Y,Z);
hold on

% Create textarrow showing the location of xglobal
annotation('textarrow',[0.44 0.50],[0.63 0.58],'TextEdgeColor',[0 0 0],...
    'TextBackgroundColor',[1 1 1],'FontSize',12,'String',{'Global Max'});
axis([-1 3.75 -3 3]);

% Plot PSO and fmincon results using a circle marker
P1 = scatter(x(1), x(2), 'o', ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor',[1 0 0],'LineWidth',1.25);

P2 = scatter(xlocal(1), xlocal(2), 's', ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor',[0 0 1],'LineWidth',1.25);

legend([P1, P2],{'PSO', 'MATLAB:fmincon()'}, 'Location','best');

title('PSO Performance with Relaxed Bounds');






