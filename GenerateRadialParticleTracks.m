function allTracks = GenerateRadialParticleTracks(type,simParams)

% parameters for exponentially distributed random track lengths
numStates = length(type);
numTracks = simParams.numTracks;
Nconst = simParams.Nconst;
meanN = simParams.meanN;
minN = simParams.minN; 
maxN = simParams.maxN;
dt = simParams.dt;
% space = simParams.space;
numSubSteps = simParams.numSubSteps;
Rindex = simParams.Rindex;

radii = [];
cumRindex = cumsum([0 Rindex]);
for i = 1:length(Rindex)
    radii = [radii; cumRindex(i) cumRindex(i+1)];
end

k = 1;
X = {}; deltaX = {};
for i = 1:numStates
    
    % particle track parameters
    D = type(i).D;
    v = 0;
    sigma = type(i).sigma;
    P = type(i).P;
          
    switch type(i).name

        case 'Normal'
            simulation = 'stateTracks = SimulateDrivenDiffusion(M,numSim,D,v,sigma,dt,numSubSteps);';
                
        case 'Normal-Flow'
            v = type(i).v;
            simulation = 'stateTracks = SimulateDrivenDiffusion(M,numSim,D,v,sigma,dt,numSubSteps);';
        
        case 'Normal-Confined'
            L = type(i).L/2;
            simulation = 'stateTracks = SimulateConfinedDiffusion(M,numSim,D,L,startpos,sigma,dt,numSubSteps);';
            
        case 'fBM'
            alpha = type(i).alpha;
            simulation = 'stateTracks = SimulatefBM(M,numSim,D,alpha,sigma,dt,numSubSteps);';
            
        otherwise
            disp(['Diffusion modes match: ' type(i).name]);
    end
        
    % simulate and pase tracks
    if Nconst == 1
        M = meanN;
        numSim = round(numTracks*P);
        if strcmp(type(i).name,'Normal-Confined')
            startpos = rand(2,numSim)*2*L - L;
        end
        eval(simulation);
        x = stateTracks.observedPositionsX;
        y = stateTracks.observedPositionsY;
        x2 = stateTracks.truePositionsX;
        y2 = stateTracks.truePositionsY;
        
        % location of particle tracks
        randomRadian = rand(numSim,1)*2*pi;
        randomLocation = rand(numSim,1)*diff(radii(i,:)) + radii(i,1) + 0;
        x0 = randomLocation.*cos(randomRadian);
        y0 = randomLocation.*sin(randomRadian);


        for j = 1:numSim
            % if flow exists, then rotate particle track to create radial flow
            if v ~= 0
                theta = randomRadian(j);
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

                % noisy tracks
                pos = [x(:,j) y(:,j)]*R;
                x1 = pos(:,1);
                y1 = pos(:,2);

                % true tracks
                pos = [x2(:,j) y2(:,j)]*R;
                x3 = pos(:,1);
                y3 = pos(:,2);
            else
                x1 = x(:,j);
                y1 = y(:,j);
                x3 = x2(:,j);
                y3 = y2(:,j);
            end
        
            X{k} = [x1+x0(j) y1+y0(j)];
            deltaX{k} = diff(X{k});
            X_true{k} = [x3+x0(j) y3+y0(j)];
            deltaX_true{k} = diff(X_true{k});
            k = k + 1;
        end
    else
        
        population = round(numTracks*P);
        N = round(RandomTrackLength(population,meanN,minN,maxN));
        N = sort(N);
        minN = N(1);
        maxN = N(end);
        for M = minN:maxN
            index = find(N==M);
            if ~isempty(index)
                numSim = length(index);
                eval(simulation);
                x = stateTracks.observedPositionsX;
                y = stateTracks.observedPositionsY;
                x2 = stateTracks.truePositionsX;
                y2 = stateTracks.truePositionsY;
                
                % location of particle tracks
                randomRadian = rand(numSim,1)*2*pi;
                randomLocation = rand(numSim,1)*diff(radii(i,:)) + radii(i,1) + space;
                x0 = randomLocation.*cos(randomRadian);
                y0 = randomLocation.*sin(randomRadian);

                for j = 1:numSim
                    % if flow exists, then rotate particle track to create radial flow
                    if v ~= 0
                        theta = randomRadian(j);
                        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

                        % noisy tracks
                        pos = [x(:,j) y(:,j)]*R;
                        x1 = pos(:,1);
                        y1 = pos(:,2);

                        % true tracks
                        pos = [x2(:,j) y2(:,j)]*R;
                        x3 = pos(:,1);
                        y3 = pos(:,2);
                    else
                        x1 = x(:,j);
                        y1 = y(:,j);
                        x3 = x2(:,j);
                        y3 = y2(:,j);
                    end

                    X{k} = [x1+x0(j) y1+y0(j)];
                    deltaX{k} = diff(X{k});
                    X_true{k} = [x3+x0(j) y3+y0(j)];
                    deltaX_true{k} = diff(X_true{k});
                    disp(num2str(k));
                    k = k + 1;
                end
            end
        end
    end
end

% store tracks
allTracks.X = X;
allTracks.deltaX = deltaX;
allTracks.X_true = X_true;
allTracks.deltaX_true = deltaX_true;
