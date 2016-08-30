function results = GenerateParticleTracks(Nindex,Dindex,Sindex,A,params)

% parameters for exponentially distributed random track lengths
dt = params.dt;
numSubSteps = params.numSubSteps; 
varyT = params.varyT;
meanT = params.meanT;
minT = params.minT;
maxT = params.maxT;

k = 1;
X = {}; X_true = {}; deltaX = {}; deltaX_true = {}; stateseq = {};
for initialState = 1:length(Dindex)

    % number of protein tracks that start in this state
    numSim = Nindex(initialState);

    % get random track length
    if varyT == 1
        T = round(RandomTrackLength(numSim,meanT,minT,maxT));
    else
        T = ones(numSim,1)*meanT;
    end
    
    for j = 1:numSim

        results = SimulateDiffusionTransitionsSingle(initialState,T(j),Dindex,Sindex,A,dt,numSubSteps);
        
        % parse tracks
        stateseq{k} = results.stateseq;
        
        X{k} = [results.observedPositionsX results.observedPositionsY];
        deltaX{k} = diff(X{k});

        X_true{k} = [results.truePositionsX results.truePositionsY];
        deltaX_true{k} = diff(X_true{k});

        k = k + 1;
    end
    
end

results.stateseq = stateseq;
results.X = X;
results.deltaX = deltaX;
results.X2 = X_true;
results.deltaX_true = deltaX_true;

