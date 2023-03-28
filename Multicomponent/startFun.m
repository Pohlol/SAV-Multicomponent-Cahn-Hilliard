function [val] = startFun(nodes,k,N)
    %possible starting function, but gets complicated very quickly
    %only works for one case :)
    if N==3
        if k==1
            val = 1/2*cos(pi/5*nodes(:,1)).*cos(pi/5*nodes(:,2))+1/2;
        elseif k==2
            val =1/2-1/2*cos(pi/5*nodes(:,1)).*cos(pi/5*nodes(:,2));
        elseif k==3
            val = zeros(length(nodes),1);
        else 
            error("There are only dimensions");
        end
    elseif N==4
        error("No starting values yet");
    elseif N==5
        error("No starting values yet");
    else
        error("dimension too high")
    end
end

