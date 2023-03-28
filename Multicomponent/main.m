%this program solves the multicomponent Cahn-Hilliard equations with 
% homogeneous Neumann boundary conditions using
%Finite Element Methode for space and the SAV Method for time


%choose the number of components
N=3;
%set the size of the time interval
T=5;
%set value for tau and epsilon
tau = 0.003;
epsilon = 0.1;
% get the number of timesteps
number_of_timesteps = floor(T/tau);

%get the nodes and elements
nodes = load("Nodes/Square_nodes5.txt");
elements = load("Nodes/Square_elements5.txt");

number_of_nodes = length(nodes);
number_of_elements = length(elements);

%get the mass and stiffness matrices
[A,M] = assembly_bulk(nodes,elements);

%precallocate space for the variables
phi = zeros(number_of_nodes,number_of_timesteps,N);
mu = zeros(number_of_nodes,number_of_timesteps,N);
r = zeros(1,number_of_timesteps);

%set the starting values for phi
%NICHT FERTIG
% for k=1:N
%     phi(:,1,k) = startFun(nodes,k,N);
% end

%random starting values
rng(5);
phi(:,1,1) = 1/2*rand(length(nodes),1);
rng(4);
phi(:,1,2) = 1/2*rand(length(nodes),1);
phi(:,1,3) = 1-phi(:,1,1)-phi(:,1,2);

%calculate the starting value for r
r(1) = sqrt(numInt(elements,nodes,number_of_elements,FF(phi(:,1,:),N)));

%calculate the solution for all timesteps
for i = 2:number_of_timesteps
    
    E_h = sqrt(numInt(elements,nodes,number_of_elements,FF(phi(:,i-1,:),N)));
    DF = DFF(phi(:,i-1,:),N);
    bigM = kron(eye(N-1),M);
    bigA = kron(eye(N-1),A);
    
    %calculates the integral on the r.h.s. in the third equation
    rhs_int = 0;
    for n1=1:N
        rhs_int = rhs_int + transpose(DF(:,n1))*M*phi(:,i-1,n1);
    end
    
    %r.h.s. for the first equation
    rhs_help = zeros((N-1)*number_of_nodes,1);
    for n2=1:N-1
         rhs_help((n2-1)*number_of_nodes+1:n2*number_of_nodes) = M*phi(:,i-1,n2);
    end

    %put together all the parts for r.h.s.
    rhs = [rhs_help;
           zeros((N-1)*number_of_nodes,1);
           tau*r(i-1)+tau/(2*E_h)*transpose(M*DF(:,N))*ones(number_of_nodes,1)-tau/(2*E_h)*rhs_int];
    
    %calculates the integral on the l.h.s. in the third equation
    lhs_int = M*DF(:,1)-M*DF(:,N);
    if N>2
        for n3 = 2:N-1
        lhs_int = cat(1,lhs_int,M*DF(:,n2)-M*DF(:,N));
        end
    end

    %projection onto H_0
    ProjDF = projecton(DF,N);

    %integral on the l.h.s. in the second equation 
    lhs_help = zeros((N-1)*number_of_nodes,1);
    for n4=1:N-1
        lhs_help((n4-1)*number_of_nodes+1:n4*number_of_nodes) = M*ProjDF(:,n4);
    end

    %put together the l.h.s.
    lhs = [bigM,tau*bigA,zeros((N-1)*number_of_nodes,1);
           -tau*epsilon*bigA,tau*bigM,-tau/(epsilon*E_h)*lhs_help;
           -tau/(2*E_h)*transpose(lhs_int),zeros(1,(N-1)*number_of_nodes),tau];

    sol = linsolve(lhs,rhs);
    for n5 = 1:N-1
        phi(:,i,n5) = sol((n5-1)*number_of_nodes+1:n5*number_of_nodes);
        mu(:,i,n5) = sol((n5+N-2)*number_of_nodes+1:(n5+N-1)*number_of_nodes);
    end
    r(i) = sol(end);
    phi(:,i,N) = ones(number_of_nodes,1) - sum(phi(:,i,1:N-1),3);
    mu(:,i,N) = -sum(mu(:,i,1:N-1),3);
    i
end

    %calculate the mass and energy at each step
    mass = zeros(N,number_of_timesteps);
    FE = zeros(1,number_of_timesteps);
    modFE = zeros(1,number_of_timesteps);
    diff = zeros(1,number_of_timesteps);
    gradPhi = zeros(1,number_of_timesteps);
for l=1:number_of_timesteps
    for n6=1:N
        mass(n6,l) = numInt(elements,nodes,number_of_elements,phi(:,l,n6));
        FE(l) = FE(l) + epsilon/2*transpose(phi(:,l,n6))*A*phi(:,l,n6);
    end
    gradPhi(l) = FE(l);
    modFE(l) = FE(l) + r(l)^2/epsilon;
    FE(l) = FE(l) + 1/epsilon*numInt(elements,nodes,number_of_elements,FF(phi(:,l,:),N));
    diff(l) = r(l)^2-numInt(elements,nodes,number_of_elements,FF(phi(:,l,:),N));

end


timesteps = zeros(1,number_of_timesteps);
    for i=1:number_of_timesteps
    timesteps(i) = (i-1)*tau;
    end
%plots the mass for each component

figure    
    for n7=1:N
    subplot(3,1,n7)
    plot(timesteps, mass(n7,:))
    h3=title(['Masse f\"ur die $',num2str(n7),'$-te Komponente: $\int_\Omega \varphi_',num2str(n7),' dx$']);
    set(h3,'Interpreter','latex','FontSize',16);
    axis([0 T mass(n7,1)-0.1 mass(n7,1)+0.1])
    end


%plots the energy and modified energy
figure
    subplot(2,1,1)
    plot(timesteps,FE)
    h1=title('Energie: $\mathcal{E}(\varphi) = \int_\Omega \frac{\varepsilon}{2}|\nabla\varphi|^2 + \frac{1}{\varepsilon}F(\varphi) dx$');
    set(h1,'Interpreter','latex','FontSize',16);
    axis([-0.1 T FE(end)-5 FE(1)+5])
    hold on
    subplot(2,1,2)
    plot(timesteps, modFE)
    h2=title('SAV-Energie: $\tilde{\mathcal{E}}(\varphi) = \int_\Omega \frac{\varepsilon}{2}|\nabla\varphi|^2 dx + \frac{1}{\varepsilon}|r|^2$');
    set(h2,'Interpreter','latex','FontSize',16);
    axis([-0.1 T modFE(end)-5 modFE(1)+5])
%     subplot(4,1,3)
%     plot(timesteps,gradPhi)
%     h3=title('Nur Gradient: $\int_\Omega \frac{\varepsilon}{2}|\varphi|^2 dx$');
%     set(h3,'Interpreter','latex','FontSize',16);
%     subplot(4,1,4)
%     plot(timesteps,diff)
%     h4=title('Differenz: $|r|^2 - \int_\Omega F(\varphi) dx$');
%     set(h4,'Interpreter','latex','FontSize',16);

%plots the evolution of phi at the chosen timepoints set in times array
lsg = -phi(:,:,3)+phi(:,:,2)+2*phi(:,:,1);

times = [0 1 2 3 4 5 6 6.25 6.3875];

    figure
    
    for k = 1 : 9
        subplot(3,3,k)
        
        time = times(k);
        time_step = floor(time/tau)+1;
        trisurf(elements, nodes(:,1), nodes(:,2),0*nodes(:,2),lsg(:,time_step));
        view([0 0 1]);
        shading interp
        h1 = title(['t = ',num2str(time)]);
        set(h1,'Interpreter','latex','FontSize',16);
        axis('equal');
        hold on 
        colorbar
        colormap("jet")
        pbaspect([1,1,1])
    end


    %uncomment to make a video of the evolution of phi
    %capture_video(number_of_timesteps,nodes,elements,phi,tau);

