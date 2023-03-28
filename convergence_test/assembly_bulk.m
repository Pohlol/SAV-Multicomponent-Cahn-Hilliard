function [A,M]=assembly_bulk(nodes,elements)

    % initialization
    N=length(nodes);
    M = zeros(N,N);
    A = zeros(N,N);

    % the basis functions on the reference element
    % x     - (1,0)
    % 1-x-y - (0,0)
    % y     - (0,1)
   Mtx_m = [1/12    1/24    1/24
            1/24    1/12    1/24
            1/24    1/24    1/12];
        
   % \pa_x v^2
   Mtxx = [1/2   -1/2     0
           -1/2    1/2     0
            0      0      0];
        
   % \pa_y v^2
   Mtxy = [0     0       0
            0    1/2    -1/2
            0   -1/2     1/2 ];

        
   % \pa_x v \pa_y v
   Mtxxy = [0    -1/2    1/2 
            0    1/2     -1/2
            0     0       0  ];


    number_of_elements=size(elements,1);

    % assembly
    for i=1:number_of_elements
        AA=[nodes(elements(i,1),1) nodes(elements(i,1),2)];
        BB=[nodes(elements(i,2),1) nodes(elements(i,2),2)];
        CC=[nodes(elements(i,3),1) nodes(elements(i,3),2)];
        X0=BB(1);
        Y0=BB(2);

        X2=AA(1);
        Y2=AA(2);

        X1=CC(1);
        Y1=CC(2);

        a11=X2-X0;
        a12=X1-X0;

        a21=Y2-Y0;
        a22=Y1-Y0;

        M_A=[a11 a12;a21 a22];

        C=inv([a11 a12; a21 a22]);
        determ=det(M_A);
        
        %local mass and stiffness matrix
        A_loc=(Mtxx*(C(1,1)^2+C(1,2)^2) + ...
              (Mtxxy+Mtxxy')*(C(1,1)*C(2,1)+C(1,2)*C(2,2)) +...
               Mtxy*(C(2,1)^2+C(2,2)^2))*abs(determ);



        M_loc=Mtx_m*abs(determ);

        % mass matrix
        M(elements(i,:),elements(i,:))=M(elements(i,:),elements(i,:))+M_loc;

        % stiffness matrix
        A(elements(i,:),elements(i,:))=A(elements(i,:),elements(i,:))+A_loc;

    end
end