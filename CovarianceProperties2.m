function [vacf C] = CovarianceProperties2(deltaX,T,trackparams)

D = trackparams.D;
N = trackparams.N;

% calculate MSD and correlation of data for all particle tracks
numTracks = length(deltaX);
dim = size(deltaX{1},2);
vacf = zeros(numTracks,D,dim);
for n = 1:numTracks
    for j = 1:dim
        delta = deltaX{n}(:,j);
        xbar = mean(delta);
%         C = (delta - xbar)*(delta - xbar)';
        C = (delta)*(delta)';
%         vacf(n,:) = vacf(n,:) + CalculateVACF(C,D)/dim;
        vacf(n,:,j) = CalculateVACF(C,D);
    end
end 

% C = cell(N,1);
% for n = 1:N
%     C{n} = zeros(D);
%     for j = 1:dim
%         C{n} = C{n} + deltaX{n}(:,j)*deltaX{n}(:,j)';
%     end
%     C{n} = C{n}/dim;
% end
vacf = vacf(:,1:T,:);
C = cell(N,1);
for n = 1:N
     C{n} = toeplitz([vacf(n,1:T) zeros(1,D-T)]);
end

