function f1 = convplots_figMCSpace

% Creates convergence plots with h on the x-axis

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;

%%

T=1;    % final time

%only works if both vectors are the same size
n_vect = 2:7;
tau_vect =.2*2.^(-2:-1:-7);


% legend location
location_text='SouthEast';%'NorthWest';%
% axis limit and reference line height
ylimits=[10^-3 10^1000];

h=[0.2 0.15 0.1 0.07 0.05 0.04 0.03 0.025 0.02];
h_vect = h(n_vect);
%% reading errors
for i=1:length(n_vect)
    meshno=n_vect(i);
  
    
    for j=1:length(tau_vect)
        tau=tau_vect(j);
                err=load(['Errors/errors_MC_n',num2str(meshno),'_tau',num2str(tau),'.txt']);
                % L^\infty errors on [0,T]
                error_phi_1(i,j)=max(err(2,:));
                error_phi_2(i,j)=max(err(3,:));
                error_mu_1(i,j)=max(err(4,:));
                error_mu_2(i,j)=max(err(5,:));
                error_r(i,j)= max(err(6,:));
          Tau_legend{i}=strcat(['$\tau=',num2str(tau_vect(i)),'$']);
    end
end

Tau_legend{i+1} = '$\mathcal{O}(h^2)$';

%% plot
f1 = figure
% x 0
% 0 0
subplot(1,5,1)
loglog_conv(h_vect,error_phi_1)
hold on;
plot(h_vect,(h_vect*1.5).^2,'-- black','LineWidth',1);
hold off; 
h1=title('$\|\varphi_1-(\varphi_1^n)\|_{L^\infty(0,T;L^2(\Omega))}$');
set(h1,'Interpreter','latex','FontSize',16);
h1=xlabel('Gitterweite ($h$)');
set(h1,'Interpreter','latex','FontSize',14);
h2=legend(Tau_legend,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.1)*min(h_vect) 1.1*max(h_vect)])
ylim(ylimits)


subplot(1,5,2)
loglog_conv(h_vect,error_phi_2)
hold on;
plot(h_vect,(h_vect*1.5).^2,'-- black','LineWidth',1);
hold off; 
h1=title('$\|\varphi_2-(\varphi_2^n)\|_{L^\infty(0,T;L^2(\Omega))}$');
set(h1,'Interpreter','latex','FontSize',16);
h1=xlabel('Gitterweite ($h$)');
set(h1,'Interpreter','latex','FontSize',14);
h2=legend(Tau_legend,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.1)*min(h_vect) 1.1*max(h_vect)])
ylim(ylimits)
% 0 x
% 0 0
subplot(1,5,3)
loglog_conv(h_vect,error_mu_1)
hold on;
plot(h_vect,(h_vect*1.9).^2,'-- black','LineWidth',1);
hold off; 
h1=title('$\|\mu_1-(\mu_1^n)\|_{L^\infty(0,T;L^2(\Omega))}$');
set(h1,'Interpreter','latex','FontSize',16);
h1=xlabel('Gitterweite ($h$)');
set(h1,'Interpreter','latex','FontSize',14);
h2=legend(Tau_legend,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.1)*min(h_vect) 1.1*max(h_vect)])
ylim(ylimits)

subplot(1,5,4)
loglog_conv(h_vect,error_mu_2)
hold on;
plot(h_vect,(h_vect*1.9).^2,'-- black','LineWidth',1);
hold off; 
h1=title('$\|\mu_2-(\mu_2^n)\|_{L^\infty(0,T;L^2(\Omega))}$');
set(h1,'Interpreter','latex','FontSize',16);
h1=xlabel('Gitterweite ($h$)');
set(h1,'Interpreter','latex','FontSize',14);
h2=legend(Tau_legend,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.1)*min(h_vect) 1.1*max(h_vect)])
ylim(ylimits)


subplot(1,5,5)
loglog_conv(h_vect,error_r)
hold on;
plot(h_vect,(h_vect*1.2).^2,'-- black','LineWidth',1);
hold off; 
h1=title('$\|r-(r^n)\|_{L^\infty(0,T)}$');
set(h1,'Interpreter','latex','FontSize',16);
h1=xlabel('Gitterweite ($h$)');
set(h1,'Interpreter','latex','FontSize',14);
h2=legend(Tau_legend,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.1)*min(h_vect) 1.1*max(h_vect)])
ylim(ylimits)




% % fig and eps
 %saveas(f1,['MC_time_Linfty.fig']);
% saveas(f1,['figures/convplot_',export_txt,'_time_Linfty.eps'], 'epsc2');



end
    

function loglog_conv(vect,err)

[n1 n2]=size(err);

symbols='sox.d+^*v><';
ms=[6 6 6 8 6 6 6 6 6 6 6];
gr=(linspace(.66,0,n1))';
colors=[gr gr gr];

for jj=1:n2
    loglog(vect,err(:,jj), ...
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
