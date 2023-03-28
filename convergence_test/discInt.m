function [intOutput] = discInt(nodes,elements, number_of_elements, fun)
%integrates a function over the domain by approximating the integral on
%every triangle
intOutput = 0;
for i=1:number_of_elements
        AA=[nodes(elements(i,1),1) nodes(elements(i,1),2)];
        BB=[nodes(elements(i,2),1) nodes(elements(i,2),2)];
        CC=[nodes(elements(i,3),1) nodes(elements(i,3),2)];
        X = [AA(1),BB(1),CC(1)];
        Y = [AA(2),BB(2),CC(2)];
        intOutput = intOutput + polyarea(X,Y)*(fun(elements(i,1))+fun(elements(i,2))+fun(elements(i,3)))/3;

end

