function local_error(nodes,elements,A,M,T,meshno,tau)
%calculates the error for phi, mu and r using a quadratic potential. 
%Therefore we can calculate an exact solution by hand which gets saved in
%the variable phi,mu and r.
%the numerical values get saved in numPhi,numMu and numR
%Then calculates the L^2 norm of the difference at the end and saves it in
%a textfile

number_of_timesteps=floor(T/tau);
number_of_nodes = length(nodes);
number_of_elements = length(elements);

%preallocate space
phi = zeros(number_of_nodes,number_of_timesteps);
mu = zeros(number_of_nodes,number_of_timesteps);
r = zeros(1,number_of_timesteps);
numPhi = zeros(number_of_nodes,number_of_timesteps);
numMu = zeros(number_of_nodes,number_of_timesteps);
numR = zeros(1,number_of_timesteps);
E_h = zeros(1,number_of_timesteps);
errors = zeros(4,number_of_timesteps);
%get the time t^n at for a timestep
time = zeros(1,number_of_timesteps);
for i=1:number_of_timesteps
time(i) = tau*(i-1);
end
%set exact values for phi,mu and r
for i=1:number_of_timesteps
    phi(:,i) = solution(time(i),nodes);
    mu(:,i) = (2*(pi/5)^2+1)*solution(time(i),nodes);
    r(i) = rExact(time(i));
end
numPhi(:,1) = solution(0,nodes);
numR(1) = sqrt(discInt(nodes,elements,number_of_elements,doubleWell(numPhi(:,1))));
for i=2:number_of_timesteps
    %calculate E_h;
    E_h(i) = discInt(nodes,elements,number_of_elements,doubleWell(numPhi(:,i-1)));
    %build the r.h.s.

%     
    rhs = [M*numPhi(:,i-1)+M*f1(nodes,time(i));
        tau*M*f2(nodes,time(i));
        tau*numR(i-1)-tau/(2*sqrt(E_h(i)))*transpose(diffDoubleWell(numPhi(:,i-1)))*M*numPhi(:,i-1)];
    %build the l.h.s.
    lhs = [M,tau*A,zeros(number_of_nodes,1);
        -tau*A, tau*M, -tau/(sqrt(E_h(i)))*M*diffDoubleWell(numPhi(:,i-1));
        transpose(-tau/(2*sqrt(E_h(i)))*M*diffDoubleWell(numPhi(:,i-1))), zeros(1,number_of_nodes),tau];
    

    sol = linsolve(lhs,rhs);
    numPhi(:,i) = sol(1:number_of_nodes);
    numMu(:,i) = sol(number_of_nodes+1:2*number_of_nodes);
    numR(i) = sol(2*number_of_nodes+1);
    
%calculate the norms


errors(1,:) = time(1:number_of_timesteps);
errors(2,:) = calcNorm(phi-numPhi,M,number_of_timesteps);
errors(3,:) = calcNorm(mu-numMu,M,number_of_timesteps);
errors(4,:) = abs(r-numR);

errors = errors(:,2:number_of_timesteps);
save(strcat(['Errors/errors_SAV_n',num2str(meshno),'_tau',num2str(tau),'.txt']), 'errors', '-ASCII');



