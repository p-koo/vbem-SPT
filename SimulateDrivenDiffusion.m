function results = SimulateDrivenDiffusion(N,numSim,D,v,sigma,dt,numSubSteps)

% simulate x and y displacements with microsteps time steps
if nargin == 6
    numSubSteps = 50;
end
microtime = dt/numSubSteps;

% generate Weiner displacements
if length(v) > 1
    displacementsX = randn(N*numSubSteps,numSim)*sqrt(2*D*microtime)+v(1)*microtime;
    displacementsY = randn(N*numSubSteps,numSim)*sqrt(2*D*microtime)+v(2)*microtime;
else
    displacementsX = randn(N*numSubSteps,numSim)*sqrt(2*D*microtime)+v*microtime;
    displacementsY = randn(N*numSubSteps,numSim)*sqrt(2*D*microtime);
end

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





