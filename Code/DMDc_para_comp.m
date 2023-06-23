function [] = DMDc_para_comp(A_tilde, B_tilde, U_x , rank, jumps ,steps,cycles )

%%% This script computes the DMDc for the 5 cycles version 
%%% for comparison with other parametric system



%hemodynamic parameter
nu = 0.004/1055;
%snap vector length
snap_vec = 5136636;
%step size in time
dt = 0.02;
% .bin file for comparison data
snapdata = '5cn_all.bin';
%number of snaps given by cycles & jumps
%snaps = floor(spc*cycles/jumps);
%matrix for input data
Y_para = zeros(4,steps);

%alternative input data (mass flow rate)
in_data_para = readmatrix('STAR_massFlowRate_schindler_5zyklen_new.txt');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read alternative FEM data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load all snapshots
id=fopen(snapdata);
X_comp = fread(id,[snap_vec 250],'double');
fclose(id);
%reduce data to the necessary snaps 
X_comp = X_comp(:,1:jumps:steps);
%how many steps for comparison in total given interval
compsteps = size(X_comp,2);

%% load relevant inflow data
for i=1:steps
        %find current inflow value
        [ ~, idx]= min( abs( in_data_para(:,1) - dt*i ) );
        m_tilde = in_data_para(idx,2);
        
        %time deriv. inflow.
        m_der = (m_tilde - in_data_para(idx-1,2))/(in_data_para(idx,1) - in_data_para(idx-1,1));
        
        %input snapshot at time t
        y = [ m_der ; m_tilde^2 ; m_tilde ; nu*m_tilde ];
        Y_para(:,i) = y ;
end

%reduced Y for DMD according to step shifts
Y_para = Y_para(:,1:jumps:steps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solution comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
x_dis = zeros(rank,compsteps);
X_dis = zeros(snap_vec,compsteps);
n_comp = zeros(compsteps,1);
%totaler Fehler
nn = 0;


%reconstruct first snapshot
x_dis(:,1)=eye(rank)/A_tilde * (U_x' * X_comp(:,2) - B_tilde*Y_para(:,1));

X_dis(:,1) = U_x * x_dis(:,1);
n_comp(1) = norm( X_comp(:,1)-X_dis(:,1) )/norm(X_comp(:,1));


%rescale ROM x to full size and compute relative error
for k = 2:compsteps
    x_dis(:,k) = A_tilde * x_dis(:,k-1) + B_tilde * Y_para(:,k-1);
    X_dis(:,k) = U_x * x_dis(:,k);
    n_comp(k) = norm( X_comp(:,k)-X_dis(:,k) )/norm(X_comp(:,k));
    nn = nn + n_comp(k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% error plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot limit of input
[ ~, idx]= min( abs( in_data_para(:,1) - steps*dt ) );

%rel. error plots
figure()
hold on
yyaxis left
plot(in_data_para(1:idx,1),in_data_para(1:idx,2))
ylabel('MassFlowRate Input  $\dot{m}$', 'interpreter','latex', 'fontsize',16)
hold off
yyaxis right
hold on
semilogy([dt+dt*jumps:dt*jumps:steps*dt],n_comp(2:end))
semilogy([0 5],[1e-1 1e-1])
semilogy([0 5],[3e-1 3e-1], 'r')
hold off
set(gca,'yscale','log')
ylabel(sprintf('rel. error for $\\Delta t=%.2f$ , rank=%i, cycles=%.1f',dt*jumps,rank, cycles), 'interpreter','latex', 'fontsize',16)
xlabel('time $t$ in $s$', 'interpreter','latex', 'fontsize',18)
title('DMDc comparison with alt. setup', 'fontsize',24)
grid minor
h=gcf;
h.WindowState = 'maximized';
path = fullfile('D:\Eigene Dokumente\Uni\BA',sprintf('para_comp_%iaus%.1f_r%i.fig',jumps,cycles,rank));
savefig(path)
close figure 1



end
