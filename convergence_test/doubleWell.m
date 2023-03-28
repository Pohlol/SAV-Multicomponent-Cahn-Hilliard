function [F] = doubleWell(phi)
%F = 1/4*(1-phi.^2).^2 + 0.0001;
%F=ones(size(phi));
F = 1/2*phi.^2 + 0.0001;
end

