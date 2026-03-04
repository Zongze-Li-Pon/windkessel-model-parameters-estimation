%% 
clc, clear, close all;

%% Section A: necessary information for estimating the parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: these information are required for solving the ODE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: flowrate of AAo (in mL/s).
data = readmatrix('flowrate.csv'); % This can be replaced by the csv data
                                   % you have.
T = data(:,1)';
Q = data(:,2)';
% Flow distribution, QPercent is the flow percentage for each branch.
QPercent = 0.2; % You may want to change this to have the flowrate for 
                % different branch. This percentage will multiply by 
                % the AAO flowrate.
Q = Q * QPercent;

% Objective function: desired max and min pressure in resultant waveform.
PMax = 50; PMin = 10; % Change these if you have different goal.

% Range of compliance for all branches. This is required when MATLAB
% searching the optimal set of Windkessel model parameters. Normally, we
% don't need to change these since these cover most of the compliance
% reported in the field.
CMax = 3; CMin = 0;

% Initial pressure for ODE solver. 
P0 = 10;

% How many periods will run for ODE. The resultant waveform may have
% transient periods according to the initial pressure you set above. If it
% is too off from the final pressure waveform, the transient periods may be
% longer. 20 periods are a good duration to eliminate all the transient
% periods. You can also check if it reach a steady period from the output
% figure.
nT = 20;

% Time step for ODE solver.
period = T(end);
dt = period / 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: these information are required relavent to the algorithm
% mentioned in the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometric resistance, obtained from the steady flow simulation.
Rgeo = 0.176160537; % Change this for your specific geometry.

%% Section B: pre-process the parameters.
% Calculate the average flowrate.
avFR = 0;
for i = 2:numel(Q)
    avFR = avFR+(Q(i)+Q(i-1))*(T(i)-T(i-1))/2;
end
avFR = avFR/T(end);
fprintf("Av. Q is: %f.\n", avFR);

% Calculate the mean pressure.
PMean = 1/3*PMax+2/3*PMin;
fprintf("Mean P is: %f.\n", PMean);

% Calculate the total resistance.
RTotal = PMean/avFR;
fprintf("Total resistance is: %f\n", RTotal);
fprintf("Geometric resistance is: %f\n", Rgeo);

% Calculate the total resistance for the outlet Windkessel model.
RTotal = RTotal-Rgeo;
fprintf("Total resistance (without Rgeo) is %f\n", RTotal);

% Calculate the desired pressure for each branch outlet.
PObjMax = PMax-avFR*Rgeo;
PObjMin = PMin-avFR*Rgeo;
fprintf("Objective max and min pressure are %f and %f\n", PObjMax, PObjMin);

% Generate a flowrate profile with n periods.
flowrate = repmat(Q(1:end-1), 1, nT);
flowrate = [flowrate,Q(end)];

% Generate a time profile with n iterations
for i = 0:nT-1
    time(1+i*(numel(T)-1):numel(T)-1+i*(numel(T)-1)) = T(1:end-1)+i*T(end);
end
time = [time,nT*T(end)];

% Find the gradient of the flowrate for Windkessel model and interplate
% them with the same x.
flowrateG = gradient(flowrate, time);
x = 0:dt:time(end);
Q1 = interp1(time, flowrate, x, 'spline');
flowrateG1 = gradient(Q1, x);

%% Section C: use pattern search to find the optimal RCR parameters.
tic

% Range setup for the search.
R1 = optimvar('R1',1, "LowerBound",0,"UpperBound",RTotal);
C = optimvar('C',1, "LowerBound",CMin,"UpperBound",CMax);

% Convert the function to optimization expression by using fcn2optimexpr.
tspan = 0:dt:time(end);
myfcn = fcn2optimexpr(@RCtoODE, R1, C, RTotal, period, x, Q1, flowrateG1, ...
                      tspan, P0);

% Express the objective function as the sum of squared differences between
% the ODE solution and Pmax and Pmin target.
obj = sqrt(sum((myfcn - [PObjMin, PObjMax]).^2));

% Create an optimization problem with the objective function obj.
prob = optimproblem("Objective",obj);

% Give an initial guess r0 for the solver and call solve. Solving problem
% using lsqnonlin.
ini.R1 = RTotal/3;
ini.C = (CMax+CMin)/2;

% Patternsearch solver
options = optimoptions('patternsearch', 'UseParallel', true,'PlotFcn', ...
                       @psplotbestf, 'MaxIterations', 1000);
[rsol, objSol, exitflag1, output1] = solve(prob, ini, "Solver", ...
                                           "patternsearch", "Options", ...
                                           options)
fprintf("\tOpt R1, R2, and C are: %f, %f, and %f, error is %f\n", rsol.R1, ...
        RTotal-rsol.R1, rsol.C, objSol);

toc

%% Section D: plot the solution.
% Plot the solution with optimal parameters.
figure;
subplot(2,1,1);
opts = odeset('Reltol', 1e-7, 'AbsTol', 1e-5, 'Stats', 'off');
[t, P] = ode45(@(t,P) WK3ODE(t, P, rsol.R1, rsol.C, RTotal, x, Q1, flowrateG1), ...
               tspan, P0, opts);
plot(t,P);

% The last period result has larger errors for Pmin due to the truncation
% of flowrate at the end of period. Therefore, use the second last period
% to calcualtion min & max.
subplot(2,1,2);
index = (t >= (nT-2)*T(end) & t<(nT-1)*T(end));
plot(t(index), P(index))
xlabel('time')
ylabel('Pressure (mmHg)')
PLast=P(index);
tLast=t(index);
fprintf("max P and min P are: %f, %f, Pressure diff in the last cycle is %f \n",...
        max(P(index)), min(P(index)), PLast(end)-PLast(1));
sumO = sqrt(sum(([min(P(index)),max(P(index))] - [PObjMin, PObjMax]).^2))

%% Section E: functions used in the above code.
function dPdt = WK3ODE(x, y, R1, C, RTotal, time, flowrate, flowrateG)
    R2 = RTotal - R1;
    if (R2 < 0)
        R2 = 0;
    end
    Q = interp1(time, flowrate, x,'spline');
    QGradient = interp1(time, flowrateG, x, 'spline');
    dPdt = (R1+R2)*Q/(R2*C)+R1*QGradient-y/(R2*C);
end

% A function that computes the ODE solution using parameters R(1)= r1, 
% R(2)= r2 & R(3)= C.
function solpt = RCtoODE(R1, C, RTotal, period, time, flowrate, flowrateG, tspan, y0)
    opts = odeset('Reltol', 1e-7, 'AbsTol', 1e-5, 'Stats', 'off');
    sol = ode45(@(t,y)WK3ODE(t, y, R1, C, RTotal, time, flowrate, flowrateG), ...
                tspan, y0, opts);
    tspan = time(end)-period*2:0.001:time(end)-period;
    solpts = deval(sol,tspan);
    solpt = [min(solpts), max(solpts)];
end