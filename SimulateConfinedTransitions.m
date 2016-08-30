function [X markovStateSeq] = SimulateConfinedTransitions(state,A,simParams)

numStates = length(state);

numTracks = simParams.numTracks;
N = simParams.N;
startLocation = simParams.startLocation;
initialMarkovState = simParams.initialMarkovState;
dt = simParams.dt;
numSubSteps = simParams.numSubSteps;

locNoise = zeros(numStates,1);
for j = 1:numStates
    locNoise(j) =  state(j).locNoise;
end

microtime = dt/numSubSteps;

% map transition probabilities onto uniform distribution
transitionCum = cumsum(A,2);

X = cell(numTracks,1); markovStateSeq = cell(numTracks,1); 
for z = 1:numTracks

    % simulate transition probabilities for each protein trajectory
    markovStateChange = rand(N,1);

    % underlying markov markovState sequences
    markovState = zeros(N,1);
    markovState(1) = initialMarkovState(z);
    for t = 2:N
        markovState(t) = find(markovStateChange(t) < transitionCum(markovState(t-1),:),1);
    end

    % set up initial confinement
    switch state(markovState(1)).mode

        case 1 % normal diffusion
            confineLocation = [];

        case 2 % confined diffusion                
            L = state(markovState(1)).L + randn(1)*state(markovState(1)).Lnoise;
            confineLocation = [-L L -L L];
        otherwise

    end

    % simulate positions
    lastmarkovState = markovState(1);
    allTruePositions = zeros(N*numSubSteps,2);
    allmarkovStateSeq = zeros(N*numSubSteps,2);
    allTruePositions(1,:) = startLocation(z,:);
    for t = 1:N
        range = 1 + (t*numSubSteps-numSubSteps+1:t*numSubSteps);

        % free diffusion
        D = state(markovState(t)).D;
        switch state(markovState(t)).mode
            
            case 1 % normal diffusion
                deltaX = randn(numSubSteps,2).*sqrt(2*repmat(D,numSubSteps,2)*microtime);
                subPositions = cumsum([allTruePositions(range(1)-1,:); deltaX]);
            
            case 2 % confined diffusion                
                L = state(markovState(t)).L + randn(1)*state(markovState(t)).Lnoise;
                
                lastPosition = allTruePositions(range(1)-1,:);

                if markovState(t) ~= lastmarkovState
                    confineLocation = [confineLocation; lastPosition(1,1)-L lastPosition(1,1)+L lastPosition(1,2)-L lastPosition(1,2)+L];
                end

                subPositions = zeros(numSubSteps+1,2);
                subPositions(1,:) = lastPosition;
                for j = 2:numSubSteps + 1

                    % propose new position
                    proposedPosition = subPositions(j-1,:) + sqrt(2*D*microtime)*randn(1,2);

                    % impose reflecting boundary conditions
                    for i = 1:2
                        if proposedPosition(i) > confineLocation(end,i*2)
                            subPositions(j,i) = 2*confineLocation(end,i*2) - proposedPosition(i); 
                        elseif proposedPosition(i) < confineLocation(end,i*2-2+1)
                            subPositions(j,i) = 2*confineLocation(end,i*2-2+1) - proposedPosition(i); 
                        elseif proposedPosition(i) >= confineLocation(end,i*2-2+1) & proposedPosition(i) <= confineLocation(end,i*2)
                            subPositions(j,i) = proposedPosition(i);
                        end
                    end
                end
                
            case 3 % driven diffusion
                deltaX = randn(numSubSteps,2).*sqrt(2*repmat(D,numSubSteps,2)*microtime) + state(markovState(t)).v*microtime;
                subPositions = cumsum([allTruePositions(range(1)-1,:); deltaX]);

            otherwise
        end
        
        % set new last markov state
        lastmarkovState = markovState(t);
        
        % save true positions and state
        allTruePositions(range,:) = subPositions(2:end,:);
        allmarkovStateSeq(range) = ones(length(range),1)*lastmarkovState;
    end

    % find true positions every dt
    truePositions = allTruePositions(1:numSubSteps:N*numSubSteps,:);

    % average over true positions to get the motion blur 
    motionPositions = zeros(N,2);
    for i = 1:N
        range = i*numSubSteps-numSubSteps+1:i*numSubSteps;
        motionPositions(i,:) = mean(allTruePositions(range,:));
    end

    % simulate noise 
    noiseDisplacement = randn(N,2).*repmat(locNoise(markovState),1,2);

    % get the observed displacements
    X{z} = motionPositions + noiseDisplacement;
    markovStateSeq{z} = markovState;
end


