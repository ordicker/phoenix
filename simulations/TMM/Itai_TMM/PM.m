% Calulates the propagation matrix in a material with thickness d and eps
function [PM] = PM(d, kz)

PM  = [exp(-1i*kz.*d) zeros(size(kz)); zeros(size(kz)) exp(1i*kz.*d)];

end