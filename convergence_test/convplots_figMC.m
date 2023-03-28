function f1 = convplots_figMC

% Creates convergence plots with tau on the x-axis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;

%%

T=1;    % final time

n_vect = 1:4;
tau_vect =.2*2.^(-1:-1:-4);


% legend location
location_text='SouthEast';%'NorthWest';%
% axis limit and reference line height
ylimits=[10^-1 10^1000];

%% reading errors
for i=1:length(n_vect)
    meshno=n_vect(i);
    Nodes=load(['Nodes/Square_nodes',num2str(meshno),'.txt']);
    dof=length(Nodes);

    DOF{i}=strcat(['dof ',num2str(dof)]);

  
    
    for j=1:length(tau_vect)
        tau=tau_vect(j);
                err=load(['Errors/errors_MC_Test_n',num2str(meshno),'_tau',num2str(tau),'.txt']);
                % L^\infty errors on [0,T]
                error_phi_1(i,j)=max(err(2,:));
%                 error_phi_2(i,j)=max(err(3,:));
%                 error_mu_1(i,j)=max(err(4,:));
%                 error_mu_2(i,j)=max(err(5,:));
%                 error_r(i,j)= max(err(6,:));
          
    end
end

    DOF{i+1}='$\mathcal{O}(\tau)$';
    DOF{i+2}='$\mathcal{O}(\tau^{1/2})$';


%% plot
f1 = figure
% x 0
% 0 0
%subplot(1,5,1)
loglog_conv(tau_vect,error_phi_1)
hold on;
plot(tau_vect,tau_vect,'-- black','LineWidth',1);
hold off; 
h1=title('$\|\phi_1-(\phi_1^n)\|_{L^\infty(L^2)}$');
set(h1,'Interpreter','latex','FontSize',18);
h1=xlabel('step size ($\tau$)');
set(h1,'Interpreter','latex','FontSize',16);
h2=legend(DOF,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.2)*min(tau_vect) 1.2*max(tau_vect)])
ylim(ylimits)


subplot(1,5,2)
loglog_conv(tau_vect,error_phi_2)
hold on;
plot(tau_vect,tau_vect,'-- black','LineWidth',1);
hold off; 
h1=title('$\|\phi_2-(\phi_2^n)\|_{L^\infty(L^2)}$');
set(h1,'Interpreter','latex','FontSize',18);
h1=xlabel('step size ($\tau$)');
set(h1,'Interpreter','latex','FontSize',16);
h2=legend(DOF,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.2)*min(tau_vect) 1.2*max(tau_vect)])
ylim(ylimits)
% 0 x
% 0 0
subplot(1,5,3)
loglog_conv(tau_vect,error_mu_1)
hold on;
plot(tau_vect,tau_vect,'-- black','LineWidth',1);
hold off; 
h1=title('$\|\mu_1-(\mu_1^n)\|_{L^\infty(L^2)}$');
set(h1,'Interpreter','latex','FontSize',18);
h1=xlabel('step size ($\tau$)');
set(h1,'Interpreter','latex','FontSize',16);
h2=legend(DOF,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.2)*min(tau_vect) 1.2*max(tau_vect)])
ylim(ylimits)

subplot(1,5,4)
loglog_conv(tau_vect,error_mu_2)
hold on;
plot(tau_vect,tau_vect,'-- black','LineWidth',1);
hold off; 
h1=title('$\|\mu_2-(\mu_2^n)\|_{L^\infty(L^2)}$');
set(h1,'Interpreter','latex','FontSize',18);
h1=xlabel('step size ($\tau$)');
set(h1,'Interpreter','latex','FontSize',16);
h2=legend(DOF,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.2)*min(tau_vect) 1.2*max(tau_vect)])
ylim(ylimits)


subplot(1,5,5)
loglog_conv(tau_vect,error_r)
hold on;
plot(tau_vect,tau_vect,'-- black','LineWidth',1);
hold off; 
h1=title('$\|r-(r^n)\|_{L^\infty}$');
set(h1,'Interpreter','latex','FontSize',18);
h1=xlabel('step size ($\tau$)');
set(h1,'Interpreter','latex','FontSize',16);
h2=legend(DOF,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.2)*min(tau_vect) 1.2*max(tau_vect)])
ylim(ylimits)




% % fig and eps
% saveas(f1,['MC_time_Linfty.fig']);
% saveas(f1,['figures/convplot_',export_txt,'_time_Linfty.eps'], 'epsc2');




end
    

function loglog_conv(vect,err)

[n1 n2]=size(err);

symbols='sox.d+^*v><';
ms=[6 6 6 8 6 6 6 6 6 6 6];
gr=(linspace(.66,0,n1))';
colors=[gr gr gr];

for jj=1:n1
    loglog(vect,err(jj,:), ...
           'LineWidth',1,...
           'Marker',symbols(jj),...
           'MarkerSize',ms(jj),...
           'Color', colors(jj,:));%'Color', 'black'); %
    if jj==1
        hold on;
    end
end
hold off;
end


