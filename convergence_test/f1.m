function [val] = f1(nodes,time)
%val = (-1+2*(pi^2)/25).*exp(-time).*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2));
val = zeros(length(nodes),1);

end

