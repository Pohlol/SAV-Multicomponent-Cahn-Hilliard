function [fun] = startFun(nodes)
  %choose your starting function for phi

    fun = cos(pi/5*nodes(:,1)).*cos(pi/5*nodes(:,2));
    %fun = cos(pi*nodes(:,1)).*cos(pi/5*nodes(:,2));
   
    %fun = cos(10*pi*nodes(:,1)).*cos(12*pi*nodes(:,2));
    %rng(5);
    %fun = 0.2*randn(length(nodes),1);
    

end

