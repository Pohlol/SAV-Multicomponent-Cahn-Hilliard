function [DF] = DF(phi)
    %Derivative of the Potential
    DF = -phi + phi.^3;
end

