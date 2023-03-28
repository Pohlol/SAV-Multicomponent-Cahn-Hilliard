function [val] = f2(nodes,time)
%val = (-2*(pi^2)/25).*exp(-time).*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2));
val = zeros(length(nodes),1);
%val = (2-2*(pi^2)/25).*exp(-time).*cos(1/5*pi*nodes(:,1)).*cos(1/5*pi*nodes(:,2))-exp(-3*time).*(cos(1/5*pi*nodes(:,1))).^3.*(cos(1/5*pi*nodes(:,2))).^3;
end

