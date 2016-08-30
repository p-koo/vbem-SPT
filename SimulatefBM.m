function results = SimulatefBM(N,numSim,D,alpha,sigma,dt,numSubSteps)

% simulate x and y displacements for 1 ms time steps
if nargin == 7
    numSubSteps = 50;
end
microtime = dt/numSubSteps;

model = 2*D;
C = zeros(N*numSubSteps,N*numSubSteps);
for n = 1:N*numSubSteps
	for m = 1:N*numSubSteps
		C(n,m) = model/2*(microtime)^alpha*(abs(n-m+1)^alpha - 2*abs(n-m)^alpha + abs(n-m-1)^alpha);
	end
end
L = chol(C,'lower');

displacementsX = L*randn(N*numSubSteps,numSim);
displacementsY = L*randn(N*numSubSteps,numSim);

% get the true x displacements
allTruePositionsX = cumsum(displacementsX);
allTruePositionsY = cumsum(displacementsY);

% find true displacements for dt
truePositionsX = allTruePositionsX(1:numSubSteps:N*numSubSteps,:);
truePositionsY = allTruePositionsY(1:numSubSteps:N*numSubSteps,:);

% average over true positions to get the motion blur 
motionPositionsX = zeros(N,numSim);
motionPositionsY = zeros(N,numSim);
for i = 1:N
    range = i*numSubSteps-numSubSteps+1:i*numSubSteps;
    motionPositionsX(i,:) = mean(allTruePositionsX(range,:));
    motionPositionsY(i,:) = mean(allTruePositionsY(range,:));
end

% simulate noise 
noiseDisplacementX = randn(N,numSim)*(sigma);
noiseDisplacementY = randn(N,numSim)*(sigma);
noiseDisplacementX = noiseDisplacementX - repmat(mean(noiseDisplacementX),N,1);
noiseDisplacementY = noiseDisplacementY - repmat(mean(noiseDisplacementY),N,1);

% get the observed displacements
observedPositionsX = motionPositionsX + noiseDisplacementX;
observedPositionsY = motionPositionsY + noiseDisplacementY;

% store results as a struct
results.observedPositionsX = observedPositionsX;
results.observedPositionsY = observedPositionsY;
results.motionPositionsX = motionPositionsX;
results.motionPositionsY = motionPositionsY;
results.truePositionsX = truePositionsX;
results.truePositionsY = truePositionsY;
results.allTruePositionsX = allTruePositionsX;
results.allTruePositionsY = allTruePositionsY;







