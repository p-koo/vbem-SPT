function results = SimulateDiffusionTransitionsSingle(initialState,T,Dindex,Sindex,A,dt,numSubSteps)

% simulate x and y displacements with microsteps time steps
if nargin == 6
    numSubSteps = 50;
end
microtime = dt/numSubSteps;

% map transition probabilities onto uniform distribution
transitionCum = cumsum(A,2);

% simulate transition probabilities for each protein trajectory
stateChange = rand(T,1);

% underlying markov state sequences
state = zeros(T,1);
state(1) = initialState;
for t = 2:T
    state(t) = find(stateChange(t) < transitionCum(state(t-1),:),1);
end

% simulate stochastic increments 
displacementsX = []; displacementsY = [];
for i = 1:T
    displacementsX = [displacementsX; randn(numSubSteps,1).*sqrt(2*repmat(Dindex(state(i,:)),numSubSteps,1)*microtime)];
    displacementsY = [displacementsY; randn(numSubSteps,1).*sqrt(2*repmat(Dindex(state(i,:)),numSubSteps,1)*microtime)];
end

% calculate true displacements every microstep
allTruePositionsX = cumsum(displacementsX);
allTruePositionsY = cumsum(displacementsY);

% find true positions every dt
truePositionsX = allTruePositionsX(1:numSubSteps:T*numSubSteps,:);
truePositionsY = allTruePositionsY(1:numSubSteps:T*numSubSteps,:);

% average over true positions to get the motion blur 
motionPositionsX = zeros(T,1);
motionPositionsY = zeros(T,1);
for i = 1:T
    range = i*numSubSteps-numSubSteps+1:i*numSubSteps;
    motionPositionsX(i,:) = mean(allTruePositionsX(range,:));
    motionPositionsY(i,:) = mean(allTruePositionsY(range,:));
end

% simulate noise 
noiseDisplacementX = randn(T,1).*Sindex(state)';
noiseDisplacementY = randn(T,1).*Sindex(state)';
noiseDisplacementX = noiseDisplacementX - repmat(mean(noiseDisplacementX),T,1);
noiseDisplacementY = noiseDisplacementY - repmat(mean(noiseDisplacementY),T,1);

% get the observed displacements
observedPositionsX = motionPositionsX + noiseDisplacementX;
observedPositionsY = motionPositionsY + noiseDisplacementY;

% store results as a struct
results.stateseq = state;
results.observedPositionsX = observedPositionsX;
results.observedPositionsY = observedPositionsY;
results.motionPositionsX = motionPositionsX;
results.motionPositionsY = motionPositionsY;
results.truePositionsX = truePositionsX;
results.truePositionsY = truePositionsY;
results.allTruePositionsX = allTruePositionsX;
results.allTruePositionsY = allTruePositionsY;

