function [val] = solutionMC_phi(time,nodes,N)
    val = zeros(length(nodes),1,N);
    val(:,1) = (-exp(-(4*(pi/5)^4+2*(pi/5)^2)*time).*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2))+1)/2;
    val(:,2) = (exp(-(4*(pi/5)^4+2*(pi/5)^2)*time).*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2))+1)/2;
end

