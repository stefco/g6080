function [s] = harmonic_partial_sum(N)

s=zeros(size(N));
for i=1:size(N(:))          % For each term in N...
    for n=1:N(i)            % ...take in up to it...
        s(i) = s(i) + 1/n;  % ...and sum the reciprocals
    end
end