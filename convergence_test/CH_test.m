function CH_test
% this program test convergence of the SAV method for Cahn-Hilliard
% equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long;


T=1;    % final time

%choose the number of finer timesteps
tau_vect =.2*2.^(-1:-1:-6); 


%size of this for loop chooses the number of meshes
for ih=1:6
    
%     figure

    
    % loading current mesh from the directory 'Nodes'
    Nodes=load(['Nodes/Square_nodes',num2str(ih),'.txt']);
    Elements=load(['Nodes/Square_elements',num2str(ih),'.txt']);


    
    
    
    
        % mass and stiffness matrix assembly
        [A,M]=assembly_bulk(Nodes,Elements);

    
    

    %% time steps loop
    for jtau=1:length(tau_vect)
        tau=tau_vect(jtau)
        local_error(Nodes,Elements,A,M,T,ih,tau)

    end
end

