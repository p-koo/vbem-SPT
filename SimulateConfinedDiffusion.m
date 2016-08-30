function results = SimulateConfinedDiffusion(N,numSim,D,L,startpos,sigma,dt,numSubSteps)


% simulate x and y displacements with microsteps time steps
if nargin == 7
    numSubSteps = 32;
end
microtime = dt/numSubSteps;

% simulate x-axis
displacementsX = randn(N*numSubSteps,numSim)*sqrt(2*D*microtime);
allTruePositionsX = zeros(N*numSubSteps,numSim);
allTruePositionsX(1,:) = startpos(1,:);
for i = 1:numSim
    for j = 1:N*numSubSteps-1
        newdisp = displacementsX(j,i);
        status = 1;
        while status == 1
            proposedPosition = allTruePositionsX(j,i) + newdisp;
            if proposedPosition <= L & proposedPosition >= -L
                allTruePositionsX(j+1,i) = proposedPosition;
                status = 0;
            elseif proposedPosition > L
                newdisp = 2*L-newdisp - 2*allTruePositionsX(j,i); 
            elseif proposedPosition < -L
                newdisp = -2*L - 2*allTruePositionsX(j,i) - newdisp;
            end
        end
    end
end

% simulate y-axis
displacementsY = randn(N*numSubSteps,numSim)*sqrt(2*D*microtime);
allTruePositionsY = zeros(N*numSubSteps,numSim);
allTruePositionsY(1,:) = startpos(2,:);
for i = 1:numSim
    for j = 1:N*numSubSteps-1
        newdisp = displacementsY(j,i);
        status = 1;
        while status == 1
            proposedPosition = allTruePositionsY(j,i) + newdisp;
            if proposedPosition <= L & proposedPosition >= -L
                allTruePositionsY(j+1,i) = proposedPosition;
                status = 0;
            elseif proposedPosition > L
                newdisp = 2*L-newdisp - 2*allTruePositionsY(j,i); 
            elseif proposedPosition < -L
                newdisp = -2*L - 2*allTruePositionsY(j,i) - newdisp;
            end
        end
    end
end

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
