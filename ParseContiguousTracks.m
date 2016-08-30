function [X trackIndex] = ParseContiguousTracks(trajectories,minLength,pixels2um)

allIndex = [];
X = {};
k = 1;
numTracks = length(trajectories);
for i = 1:numTracks

    t = trajectories{i}(:,1);
    N = length(t);

    index = [0; find(diff(t)~=1); N];
    if ~isempty(index)
        for j = 1:length(index)-1
            range = index(j)+1:index(j+1);
            X{k} = trajectories{i}(range,:);
            allIndex = [allIndex; i];
            k = k + 1;
        end
    else
        X{k} = trajectories{i};
        k = k + 1;
        allIndex = [allIndex; i];
    end
end

trackIndex = [];
X2 = {};
k = 1;
for i = 1:length(X)
    if size(X{i},1) > minLength
        X2{k} = X{i};
        trackIndex = [trackIndex; allIndex(i)];
        k = k + 1;
    end
end
data = X2;

X = cell(length(data),1);
for n = 1:length(data)
    x_raw = data{n}(:,4).*data{n}(:,6)*pixels2um;
    y_raw = data{n}(:,5).*data{n}(:,7)*pixels2um;
    X{n} = [x_raw y_raw];
end



