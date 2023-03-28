function [val] = rExact(time)
%gives exact r for F = 1/2phi^2+0.0001 and startfunction cos(x)*cos(y)
%val = sqrt(25/2*exp(-2*time)+ 0.01);
%f1 and f2 = 0 
%val = sqrt(25/2*exp(-(8*(pi/5)^4+4*(pi/5)^2)*time)+ 0.01);
%für MC
val = sqrt(25/4*exp(-(8*(pi/5)^4+4*(pi/5)^2)*time)+ 25.01);
%Für F=(1-phi^2)^2
%val = sqrt(25.01-25/2*exp(-2*time)+225/64*exp(-4*time));
end

