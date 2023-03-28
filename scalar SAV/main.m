%this program calculates a numerical solution to the Cahn-Hilliard
%equations in 2D using the Finite Element Method and the SAV method


%set the size of the time interval
T = 1;

%set the size the timestep
tau = 0.01;
number_of_timesteps = floor(T/tau);

%set the interfacial thickness parameter
epsilon = 1;
%set seed for random Data
rng(5);


% this can be used to calculate a mesh with mesh size h
% fd = @(p) drectangle(p,-1,1,-1,1);
% h = 0.05
% [nodes,elements] = distmesh2d( fd, @huniform, h , [-1,-1;1,1], [-1,-1;-1,1;1,-1;1,1] );
% nodes = 5*nodes;

%or just load nodes and elements to save time
nodes = load("Nodes/Square_nodes5.txt");
elements = load("Nodes/Square_elements5.txt");
[A,M] = assembly_bulk(nodes,elements);

number_of_nodes = length(nodes);
number_of_elements = length(elements);


%preallocate space for the variables
phi = zeros(number_of_nodes,number_of_timesteps);
mu = zeros(number_of_nodes,number_of_timesteps);
r = zeros(1,number_of_timesteps);


%set the starting values for phi and r
phi(:,1) = startFun(nodes);
r(1) = sqrt(numInt(elements,nodes,number_of_elements,F(phi(:,1))));

for i=2:number_of_timesteps
    E_h = sqrt(numInt(elements,nodes,number_of_elements,F(phi(:,i-1))));

    rhs_help = transpose(M*DF(phi(:,i-1)))*phi(:,i-1);
    rhs = [M*phi(:,i-1);
            zeros(number_of_nodes,1);
            tau*r(i-1)-tau/(2*E_h)*rhs_help];
    lhs_help = 1/(2*E_h)*M*DF(phi(:,i-1));
    lhs = [M,tau*A,zeros(number_of_nodes,1);
           -epsilon*tau*A, tau*M , -tau/epsilon*2.*lhs_help;
           transpose(-tau*lhs_help),zeros(1,number_of_nodes),tau];
    

    %use the fact that you get a linear system and solve it with linsolve
    %to get the solution
    sol = linsolve(lhs,rhs);
    phi(:,i) = sol(1:number_of_nodes);
    mu(:,i) = sol(number_of_nodes+1:2*number_of_nodes);
    r(i) = sol(2*number_of_nodes+1);
    
end


FE = zeros(1,number_of_timesteps);
modFE = zeros(1,number_of_timesteps);
mass = zeros(1,number_of_timesteps);
diff = zeros(1,number_of_timesteps);
for i=1:number_of_timesteps
    %calculate the energy and mass at each step
    
    FE(i) = 1/epsilon*numInt(elements,nodes,number_of_elements,F(phi(:,i))) + epsilon/2*transpose(phi(:,i))*A*phi(:,i);
    mass(i) = numInt(elements,nodes,number_of_elements, phi(:,i));
    %for the SAV system you only have energy dissipation for modFE, but if
    %tau and h are small enough FE and modFE are very close together. see
    %diff
    modFE(i) = 1/epsilon*r(i)^2 + epsilon/2*transpose(phi(:,i))*A*phi(:,i);
    %calculate the difference between r and int F(phi)
    diff(i) = sqrt(numInt(elements,nodes,number_of_elements,F(phi(:,i)))) - r(i);
end

%plots the evolution of phi at the timepoints set in the times array
times = [0  0.05 0.1 0.2 0.3 0.4 0.5 0.7 T-tau];

    figure
    
    for k = 1 : 9
        subplot(3,3,k)
        
        time = times(k);
        time_step = floor(time/tau)+1;
        trisurf(elements, nodes(:,1), nodes(:,2),0*nodes(:,2),phi(:,time_step));
        view([0 0 1]);
        shading interp
        h1 = title(['t = ',num2str(time)]);
        set(h1,'Interpreter','latex','FontSize',16);
        axis('equal');
        pbaspect([5 5 1])
        hold on 
        colorbar
    end

    

    timesteps = zeros(1,number_of_timesteps);
    for i=1:number_of_timesteps
    timesteps(i) = (i-1)*tau;
    end
    
    %plots the energy, modified energy and mass
    figure
    subplot(3,1,1)
    plot(timesteps,FE)
    h1=title('Energie: $\mathcal{E}(\varphi) = \int_\Omega \frac{\epsilon}{2}|\nabla\varphi|^2 + \frac{1}{\epsilon}F(\varphi) dx$');
    set(h1,'Interpreter','latex','FontSize',16);
    hold on
    subplot(3,1,2)
    plot(timesteps, modFE)
    h2=title('SAV-Energie: $\tilde{\mathcal{E}}(\varphi) = \int_\Omega \frac{\epsilon}{2}|\nabla\varphi|^2 dx + \frac{1}{\epsilon}|r|^2$');
    set(h2,'Interpreter','latex','FontSize',16);
    subplot(3,1,3)
    plot(timesteps, mass)
    h3=title('Masse: $\int_\Omega \varphi dx$');
    set(h3,'Interpreter','latex','FontSize',16);
    axis([0 T mass(1)-0.5 mass(1)+0.5])
    

    %uncommment to get a video of the evolution of phi
    %capture_video(number_of_timesteps,nodes,elements,phi,tau);

