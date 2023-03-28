function [F] = F(phi)
    %Polynomial Potential for CHE plus a small constant to keep it positive
    F = 1/4.*(1-phi.^2).^2+0.0001;
end

