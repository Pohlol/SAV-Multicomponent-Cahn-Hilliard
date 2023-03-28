function [val] = numInt(elements,nodes,number_of_elements, fun)
%calculate the integral of a function numerically by using the trapezoidal
%rule on every triangle of the mesh
val = 0;
for i =1:number_of_elements
    AA=[nodes(elements(i,1),1) nodes(elements(i,1),2)];
    BB=[nodes(elements(i,2),1) nodes(elements(i,2),2)];
    CC=[nodes(elements(i,3),1) nodes(elements(i,3),2)];
    X = [AA(1) BB(1) CC(1)];
    Y = [AA(2) BB(2) CC(2)];
    val = val + polyarea(X,Y)*1/3*(fun(elements(i,1))+fun(elements(i,2)) +fun(elements(i,3)));
end
end

