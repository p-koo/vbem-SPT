function [trajectory splitX splitIndex splitState] = SplitTracks(X,splitLength,rotate,stateSeq)

if nargin == 4
    includeState = 1;
else
    includeState = 0;
end

if nargin < 3
    rotate = 0;
end
    
n = 1;
splitX = {};
splitState = {};
splitIndex = [];
trajectory = struct([]);
for i = 1:length(X)
    
    N = length(X{i});
    k = 1;
    status = 1;
    trajectory(i).xbar = []; 
    trajectory(i).vacf = [];
    while status == 1
        range = k*splitLength-splitLength+2-k:k*splitLength+1-k;
        
        if range(end) <= N

            % rotate positions to principal axis
            pos = X{i}(range,:);            
            if rotate == 1
                C = cov(pos);
                [vec val] = eig(C);
                pos2 = (vec*pos')';            
            else
                pos2 = pos;
            end
            
            % calculate properties
            deltaX = diff(pos2);
            [xbar vacf] = CovarianceProperties(deltaX);

            % store properties
            trajectory(i).X{k} = pos2;
            trajectory(i).deltaX{k} = deltaX;
            
            trajectory(i).xbar = [trajectory(i).xbar; xbar];            
            trajectory(i).vacf = [trajectory(i).vacf; vacf];
            
            splitX{n} = deltaX;
            splitIndex = [splitIndex; i];
            if includeState == 1
                splitState{n} = stateSeq{i}(range);
            end
            
            n = n + 1;
            k = k + 1;
        else
            if range(1) <= N
                trajectory(i).excess = X{i}(range(1):end);
            else
                trajectory(i).excess = [];
            end
            status = 0;
        end
    end        
end


