function local_errorMC(nodes,elements,A,M,T,meshno,tau)

%choose the number of components
N=3;



epsilon = 1;
% get the number of timesteps
number_of_timesteps = floor(T/tau);


number_of_nodes = length(nodes);
number_of_elements = length(elements);


%precallocate space for the variables
phi = zeros(number_of_nodes,number_of_timesteps,N);
mu = zeros(number_of_nodes,number_of_timesteps,N);
r = zeros(1,number_of_timesteps);

%get the time t^n at for a timestep
time = zeros(1,number_of_timesteps);
for i=1:number_of_timesteps
time(i) = tau*(i-1);
end
%set exact values for phi,mu and r
for i=1:number_of_timesteps
    phi_exact(:,i,:) = solutionMC_phi(time(i),nodes,N);
    mu_exact(:,i,:) = solutionMC_mu(time(i),nodes,N);
    %mu(:,i) = solution(time(i),nodes);
    r_exact(i) = rExact(time(i));
end


%set the starting values for phi
%phi(:,1,3) = zeros(number_of_nodes,1);
phi(:,1,1) = (-cos(pi/5*nodes(:,1)).*cos(pi/5*nodes(:,2))+1)/2;
phi(:,1,2) = (cos(pi/5*nodes(:,1)).*cos(pi/5*nodes(:,2))+1)/2;
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
    
end
% numPhi = phi(:,:,2) - phi(:,:,1);
% numMu = mu(:,:,2) -mu(:,:,1);

%calculate the norms
errors = zeros(6,number_of_timesteps);

errors(1,:) = time(1:number_of_timesteps);
errors(2,:) = calcNorm(phi_exact(:,:,1)-phi(:,:,1),M,number_of_timesteps);
errors(3,:) = calcNorm(phi_exact(:,:,2)-phi(:,:,2),M,number_of_timesteps);
errors(4,:) = calcNorm(mu_exact(:,:,1)-mu(:,:,1),M,number_of_timesteps);
errors(5,:) = calcNorm(mu_exact(:,:,2)-mu(:,:,2),M,number_of_timesteps);
errors(6,:) = abs(r_exact-r);
errors = errors(:,2:number_of_timesteps);

save(strcat(['Errors/errors_test_MC_n',num2str(meshno),'_tau',num2str(tau),'.txt']), 'errors', '-ASCII');
