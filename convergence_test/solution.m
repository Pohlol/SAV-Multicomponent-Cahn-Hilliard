function [val] = solution(time,nodes)
    
    %val = exp(-time)*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2));
    %exact solution for f_1 =0 and f_2 = 0 and F = x^2+c_0
    val = exp(-(4*(pi/5)^4+2*(pi/5)^2)*time)*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2));
end

