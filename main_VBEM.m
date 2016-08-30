clear all;
clc;
close all;

%% diffusive state parameters

numTracks = 20000;
N = 20;
% dt = .032;
dt = .032;

Dindex = [.03 .15 .45];
Sindex = [.05 .05 .05];
Lindex = [.1 0 0];
LnoiseIndex = [.05  0 0];
Aindex = [1  1 .7];
Vindex = [0  0 0];
Pindex = [.3 .4 .3];

numStates = length(Dindex);
A = InitializeTransitionMatrix(numStates,1);


%% Simulate confined diffusion transitions

dmode = ones(1,numStates);
index = find(Lindex ~= 0);
dmode(index) = 2;
index = find(Vindex ~= 0);
dmode(index) = 3;
index = find(Aindex ~= 1);
dmode(index) = 4;

state = struct;
for i = 1:numStates
    state(i).mode = dmode(i);
    state(i).D = Dindex(i);
    state(i).locNoise = Sindex(i);
    state(i).L = Lindex(i);
    state(i).Lnoise = LnoiseIndex(i);
    state(i).v = Vindex(i);
    state(i).A = Aindex(i);
end

% simulation parameters
simParams.numTracks = numTracks;
simParams.N = N;
simParams.startLocation = zeros(simParams.numTracks,1);
simParams.dt = dt;
simParams.numSubSteps = 32;
Nindex = round(simParams.numTracks*Pindex);

% initial state
initialState = [];
for i = 1:numStates
    initialState = [initialState ones(1,Nindex(i))*i];
end
if length(initialState) < simParams.numTracks
    initialState = [ones(1,simParams.numTracks-length(initialState)) initialState];
end
simParams.initialMarkovState = initialState;

% simulate particle trajectories
[X stateSeq] = SimulateConfinedTransitions(state,A,simParams);

%%

% split tracks into equally spaced bin sizes
rotate = 0;
splitLength = N;
[trajectory splitX splitIndex splitState] = SplitTracks(X,splitLength,rotate,stateSeq);

% track parameters
T = 2; %splitLength-1;
trackparams.R = 1/6;
trackparams.dt = simParams.dt;
trackparams.dim = size(splitX{1},2);
trackparams.N = length(splitX);
trackparams.D = size(splitX{1},1);
[trackparams.vacf trackparams.C] = CovarianceProperties2(splitX,T,trackparams);

%%

% set the options for VBEM
params.maxIter = 1000;
params.lbThresh = 1;

% priors
priors.alpha0 = 1;
priors.invW0 = eye(trackparams.D)*.5;
priors.v0 = trackparams.D + 100;

% initial guess
K = 20;
minD = 1e-6;
maxD = 1;
sigma = .05;
random = 0;
[Nk Sk] = InitialGuess(K,minD,maxD,sigma,random,trackparams);

% Maximization step
est = Maximization(Nk,Sk,priors);

% Useful Expectations
E = UsefulExpectations(est,splitX,trackparams);

% Expectation step
rnk = Expectation(E,trackparams);

maxL = -Inf;
for i = 1:params.maxIter
    disp(['iteration = ' num2str(i)]);

    % sufficient statistics of obeserved data set evaluated w.r.t. responsibilities
    [Nk Sk] = SufficientStatistics(rnk,trackparams);

    % Maximization step
    est = Maximization(Nk,Sk,priors);

    % Useful Expectations
    E = UsefulExpectations(est,splitX,trackparams);

    % Expectation step
    rnk = Expectation(E,trackparams);
    sum(rnk)/trackparams.N

    % Variational lower bound step
    L = VariationalLowerBound(rnk,Nk,Sk,priors,est,E,trackparams);    
    maxL = [maxL; L];
    
    disp(num2str(maxL(end)-maxL(end-1)));
    if i > 10
        if maxL(end)-maxL(end-1) < params.lbThresh
            break;
        end
    end    
end
    
%%


