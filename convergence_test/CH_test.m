function CH_test
% C-H test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
format long;


T=1;    % final time


tau_vect =.2*2.^(-1:-1:-6); 


%% 
for ih=1:6
    
%     figure

    
    % loading current mesh from the directory 'Nodes'
    Nodes=load(['Nodes/Square_nodes',num2str(ih),'.txt']);
    Elements=load(['Nodes/Square_elements',num2str(ih),'.txt']);


    
    
    
    
        %mass and stiffness matrix assembly
        [A,M]=assembly_bulk(Nodes,Elements);

    
    

    %% time steps loop
    for jtau=1:length(tau_vect)
        tau=tau_vect(jtau)
        local_errorMC(Nodes,Elements,A,M,T,ih,tau)

    end
end


