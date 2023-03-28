function [val] = solutionMC_mu(time,nodes,N)
    val = zeros(length(nodes),1,N);
    val(:,1) = (-(1+2*pi^2/25)*exp(-(4*(pi/5)^4+2*(pi/5)^2)*time).*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2)))/2 +1/6;
    val(:,2) = ((1+2*pi^2/25)*exp(-(4*(pi/5)^4+2*(pi/5)^2)*time).*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2)))/2 +1/6;
    val(:,3) = -1/3;
end

