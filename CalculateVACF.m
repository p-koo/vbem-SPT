function vacf = CalculateVACF(C,D)

subC = zeros(D);
for i = 1:D
    subC(i,1:(D-i+1)) = C(i,i:end);
end

vacf = sum(subC)./(D:-1:1);
