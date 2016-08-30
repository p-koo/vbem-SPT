function [xbar vacf] = CovarianceProperties(deltaX)

[D dim] = size(deltaX);

% calculate empirical drift for each protein track
xbar = 0;
for j = 1:dim
    xbar = xbar + mean(deltaX(:,j))/dim;
end

% calculate empirical vacf for each protein track
vacf = zeros(1,D);
for j = 1:dim
%     delta = (deltaX(:,j)-mean(deltaX(:,j)));
    delta = deltaX(:,j);
    C = delta*delta';
    vacf = vacf + CalculateVACF(C,D);
end 

