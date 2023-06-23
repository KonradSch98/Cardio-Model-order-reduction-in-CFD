function [] = DMDc_para(cycles,jumps,rank)
%%% This script computes the DMDc for the 5 cycles version 
%%% it is a stripped down version of the large DMD implementation
%%% it calles the function for the parametric projection



tic
%hemodynamic parameter
nu = 0.004/1055;
%step size in time
dt = 0.02;
%snap vector length
snap_vec = 5136636;
%steps per cycle
spc = 50;
%how many cycles you want to use
%cycles = 2.2;
%number of total steps using in dataset
steps = 250;
%how many you want to skip
%jumps = 1;
%number of snaps given by cycles & jumps
snaps = floor(spc*cycles/jumps);
%rank = snaps-1;
% .bin file for original data
snapdata = '5c_all.bin';
%matrix for input data
Y = zeros(4,steps);

%original input data (mass flow rate)
in_data = readmatrix('STAR_massFlowRate_schindler_5zyklen.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read FEM data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load all snapshots
id=fopen(snapdata);
X_comp = fread(id,[snap_vec 250],'double');
fclose(id);
%reduce data to the necessary snaps 
X_comp = X_comp(:,1:jumps:steps);
%how many steps for comparison in total given interval
compsteps = size(X_comp,2);
%select snapshots for decomposition
X = X_comp(:,1:snaps-1);
X_p = X_comp(:,2:snaps);

%% load relevant inflow data

for i=1:steps
        %find current inflow value
        [ ~, idx]= min( abs( in_data(:,1) - dt*i ) );
        m_tilde = in_data(idx,2);
        
        %time deriv. inflow.
        m_der = (m_tilde - in_data(idx-1,2))/(in_data(idx,1) - in_data(idx-1,1));
        
        %input snapshot at time t
        y = [ m_der ; m_tilde^2 ; m_tilde ; nu*m_tilde ];
        Y(:,i) = y;
end

%reduced Y for DMD according to step shifts
Y = Y(:,1:jumps:steps);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculations

%SVD of combined matrix 
[U_o,S_o,V_o] = svd( [X;Y(:,1:snaps-1)] , 'econ');

%singular value plot
figure()
ylabel('($\Sigma$)', 'interpreter','latex', 'fontsize',14)
title('singularvalues for DMD set','fontsize',17)
semilogy([1:1:snaps-1],diag(S_o))
annotation('textbox', [0.78,0.79, 0.15, 0.15], 'String', sprintf("rank = %i\njumps = %i\ncycles = %.1f",rank,jumps,cycles),'FitBoxToText','on','FontSize',10,'BackgroundColor','w')
path = fullfile('D:\Eigene Dokumente\Uni\BA\',sprintf('SW_para_%iaus%.1f_r%i.fig',jumps,cycles,rank));
savefig(path)
close figure 1

%SVD of X_p
[U_x,~,~] = svd( X_p , 'econ');

%truncation
if rank ~= snaps-1
    U_o = U_o(:,1:rank);
    U_x = U_x(:,1:rank);
    S_o = S_o(1:rank,1:rank);
    V_o= V_o(:,1:rank);
end
U_ot1 = U_o(1:end-4, : )';
U_ot2 = U_o(end-3:end, : )';

%compute approximations and their eigenvalues
A_tilde = U_x' * X_p * V_o * eye(rank)/S_o * U_ot1 * U_x ;
B_tilde = U_x' * X_p * V_o * eye(rank)/S_o * U_ot2 ;
[~,Lamb] = eig(A_tilde);

%compute the eigenvectors /  dmd_modes 'Phi'
%U_hat = U_ot1 * U_x;
%Nu = X_p * V_o * eye(rank)/S_o * U_hat * W;


clear U_o U_ot1 U_ot2 V_o X_p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solution comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
x_dis = zeros(rank,compsteps);
X_dis = zeros(snap_vec,compsteps);
n_pro = zeros(compsteps,1);
%totaler fehler
n = 0;


%reconstruct first snapshot
x_dis(:,1)=eye(rank)/A_tilde*(U_x' * X(:,2) - B_tilde*Y(:,1));

clear X

%first snapshot
X_dis(:,1) = U_x * x_dis(:,1);
n_pro(1) = norm( X_comp(:,1)-X_dis(:,1) )/norm(X_comp(:,1));


%rescale ROM x to full size and compute relative error
for i = 2:compsteps
    x_dis(:,i) = A_tilde * x_dis(:,i-1) + B_tilde * Y(:,i-1);
    X_dis(:,i) = U_x * x_dis(:,i);
    n_pro(i) = norm( X_comp(:,i)-X_dis(:,i) )/norm(X_comp(:,i));
    n = n + n_pro(i);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% error plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot limit of input
[ ~, idx]= min( abs( in_data(:,1) - steps*dt ) );

%rel. error plots
figure()
hold on
yyaxis left
plot(in_data(1:idx,1),in_data(1:idx,2))
ylabel('MassFlowRate Input  $\dot{m}$', 'interpreter','latex', 'fontsize',16)
hold off
yyaxis right
hold on
semilogy([dt+dt*jumps:dt*jumps:steps*dt],n_pro(2:end))
semilogy([0 5],[1e-1 1e-1])
semilogy([0 5],[3e-1 3e-1], 'r')
hold off
set(gca,'yscale','log')
ylabel(sprintf('rel. error for $\\Delta t=%.2f$ , rank=%i, cycles=%.1f',dt*jumps,rank, cycles), 'interpreter','latex', 'fontsize',16)
xlabel('time $t$ in $s$', 'interpreter','latex', 'fontsize',18)
title('DMDc projection', 'fontsize',24)
grid minor
h=gcf;
h.WindowState = 'maximized';
path = fullfile('D:\Eigene Dokumente\Uni\BA',sprintf('para_%iaus%.1f_r%i.fig',jumps,cycles,rank));
savefig(path)
close figure 1




%eigenvalue plots
figure()
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'-') % plot unit circle
xlabel('Real($\Lambda(A)$)', 'interpreter','latex', 'fontsize',14)
ylabel('Imag($\Lambda(A)$)', 'interpreter','latex', 'fontsize',14)
title('eigenvalues for DMDc set','fontsize',17)
hold on, grid on
scatter(real(diag(Lamb)), imag(diag(Lamb)),'o')
axis([-1.1 1.1 -1.1 1.1]);
annotation('textbox', [0.78,0.79, 0.15, 0.15], 'String', sprintf("rank = %i\njumps = %i\ncycles = %.1f",rank,jumps,cycles),'FitBoxToText','on','FontSize',10,'BackgroundColor','w')
path = fullfile('D:\Eigene Dokumente\Uni\BA',sprintf('EW_para_%iaus%.1f_r%i.fig',jumps,cycles,rank));
savefig(path)
close figure 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% call function for parametric projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DMDc_para_comp(A_tilde, B_tilde, U_x , rank, jumps ,steps,cycles )


end