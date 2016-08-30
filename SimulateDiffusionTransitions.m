function results = SimulateDiffusionTransitions(numTracks,T,A,Nindex,Dindex,Sindex,dt,numSubSteps)

% simulate x and y displacements with microsteps time steps
if nargin == 6
    numSubSteps = 50;
end
microtime = dt/numSubSteps;

% map transition probabilities onto uniform distribution
transitionCum = cumsum(A,2);

% simulate transition probabilities for each protein trajectory
stateChange = rand(T,numTracks);

% set up initial state
initialStates = []; 
for i = 1:size(A,2)
    initialStates = [initialStates ones(1,Nindex(i))*i];
end

% underlying markov state sequences
state = zeros(T,numTracks);
state(1,:) = initialStates(1:numTracks);
for n = 1:numTracks
    for t = 2:T
        state(t,n) = find(stateChange(t,n) < transitionCum(state(t-1,n),:),1);
    end
end

% simulate stochastic increments 
displacementsX = []; displacementsY = [];
for i = 1:T
    displacementsX = [displacementsX; randn(numSubSteps,numTracks).*sqrt(2*repmat(Dindex(state(i,:)),numSubSteps,1)*microtime)];
    displacementsY = [displacementsY; randn(numSubSteps,numTracks).*sqrt(2*repmat(Dindex(state(i,:)),numSubSteps,1)*microtime)];
end


% calculate true displacements every microstep
allTruePositionsX = cumsum(displacementsX);
allTruePositionsY = cumsum(displacementsY);

% find true positions every dt
truePositionsX = allTruePositionsX(1:numSubSteps:T*numSubSteps,:);
truePositionsY = allTruePositionsY(1:numSubSteps:T*numSubSteps,:);

% average over true positions to get the motion blur 
motionPositionsX = zeros(T,numTracks);
motionPositionsY = zeros(T,numTracks);
for i = 1:T
    range = i*numSubSteps-numSubSteps+1:i*numSubSteps;
    motionPositionsX(i,:) = mean(allTruePositionsX(range,:));
    motionPositionsY(i,:) = mean(allTruePositionsY(range,:));
end

% simulate noise 
noiseDisplacementX = randn(T,numTracks).*(Sindex(state))';
noiseDisplacementY = randn(T,numTracks).*(Sindex(state))';
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

