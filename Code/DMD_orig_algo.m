function [] = DMD_orig_algo(cycles,jumps,rank)
%%% This function computes the DMD for the 5 cycle version 
%%%
%%% input -   used number of cycles for calculations
%%%           jumps says every which snap is used 
%%%           --> 1 means every snap is used, 
%%%                2 means every second is used
%%%           rank gives the truncation rank - 
%%%           --> activate the "snap-1" deinition below if you dont want to
%%%               truncate
%%%
%%% calculates also the spots in between the given snapshots if FULL==1
%%% ---> limited by FVM calculation if error has to be computed
%%% works only in the discretized scheme like the given snapshots if FULL==0
%%% also plot the interpolated snapshots if FULL_Plot==1 
%%% uses the ROM version if ROM==1
%%% uses the time-domain version if ROM==0
%%%
%%% loads snap data from .bin file
%%% ---> this hase to be in the same directory as this script



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definitions

%step size in time
dt = 0.02;
%steps per cycle
spc = 50;
%snap vector length
snap_vec = 5136636;
%how many cycles you want to use
%cycles = ;
%number of total steps using from dataset
steps = 250;
%how many you want to skip
%jumps = ;
%number of snaps given by cycles & jumps
snaps = floor(spc*cycles/jumps);
%rank = snaps-1;
% .bin file for data
snapdata = '5c_all.bin';
% script parameters
FULL = 1;
FULL_Plot = 1;
ROM = 0;
%discretized time - only necessary if FULL == 1
time = [0:dt:steps/spc-dt];


% check if script parameters are consistent
if FULL == 1 && ROM == 1 
    error('Parameters are inconsistent');
end
if FULL_Plot == 1 && FULL == 0 
    error('Parameters are inconsistent');
end


tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preperations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(snap_vec,snaps-1);
X_p = zeros(snap_vec,snaps-1);
%load all snapshots
id=fopen(snapdata);
X_comp = fread(id,[snap_vec 250],'double');
fclose(id);

%reduce data to the necessary snaps 
if FULL == 1    
    X_comp = X_comp(:,1:1:steps);
    % select snapshots for decomposition
    X = X_comp(:,1:jumps:(jumps*snaps)-jumps);
    X_p = X_comp(:,1+jumps:jumps:(jumps*snaps));

elseif FULL == 0
    X_comp = X_comp(:,1:jumps:steps);
    X = X_comp(:,1:snaps-1);
    X_p = X_comp(:,2:snaps);

end

%how many steps for comparison in total given interval
compsteps = size(X_comp,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculations

[U,S,V] = svd(X, 'econ');

% singular value plot
figure()
ylabel('($\Sigma$)', 'interpreter','latex', 'fontsize',14)
title('singularvalues for DMD set','fontsize',17)
semilogy([1:1:snaps-1],diag(S))
annotation('textbox', [0.78,0.79, 0.15, 0.15], 'String', sprintf("rank = %i\njumps = %i\ncycles = %.1f",rank,jumps,cycles),'FitBoxToText','on','FontSize',10,'BackgroundColor','w')
path = fullfile('D:\Eigene Dokumente\Uni\BA\',sprintf('SW_%iaus%.1f_r%i_%i_%i_%i.fig',jumps,cycles,rank,ROM,FULL,FULL_Plot));
savefig(path)
close figure 1


%only truncate if necessary
if rank ~= snaps-1
    U = U(:,1:rank);
    S = S(1:rank,1:rank);
    V= V(:,1:rank);
end

%compute approximation A and its eigenvalues
A_tilde = U' * X_p * V * eye(rank)/S ;
[W,Lamb] = eig(A_tilde);


%% compute the eigenvectors /  dmd_modes 'Phi'
if ROM == 0 
    Nu = X_p * V * eye(rank)/S * W;
    %compute scaling value
    b = eye(rank)/(W' * W) * W' * U' * X(:,1);
    %rescale eigenvalues
    omega = log(diag(Lamb))/(dt*jumps);

    clear U V A_tilde X X_p S
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solution evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ROM == 1
    %initialize
    x_rom = zeros(rank,compsteps-1);
    X_rec = zeros(snap_vec,compsteps-1);
    err_vec = zeros(compsteps-1,1);    
    %summed error
    n = 0;

    %first snapshot
    x_rom(:,1) = U' * X_comp(:,2);
    X_rec(:,1) = X_comp(:,2);
    err_vec(1) = norm( X_comp(:,2)-X_rec(:,1) )/norm(X_comp(:,2));
    n = n + err_vec(1);
    
    %steps and relative error
    for k = 2:compsteps-1
        x_rom(:,k) = A_tilde * x_rom(:,k-1);
        X_rec(:,k) = U * x_rom(:,k);
        err_vec(k) = norm( X_comp(:,k+1)-X_rec(:,k) )/norm(X_comp(:,k+1));
        n = n + err_vec(k);
    end


elseif ROM == 0
    if FULL == 0
        %initialize
        X_rec = zeros(snap_vec,compsteps-1);
        err_vec = zeros(compsteps-1,1);
        time = [0:dt*jumps:steps/spc-dt];
        %summed error
        n = 0;
        
        for j=1:compsteps-1
            %reconstruction and projection
            dyn = b.*exp(omega*time(j));
            X_rec(:,j) = Nu*dyn;
            %rel. error
            err_vec(j) = norm( X_comp(:,j+1)-X_rec(:,j) )/norm(X_comp(:,j+1));
            n = n + err_vec(j);
        end
        
    elseif FULL == 1
        X_rec = zeros(snap_vec,compsteps-jumps);
        err_vec = zeros(compsteps-jumps,1);
        %summed error
        n = 0;
        n_full = 0;

        count=jumps-1;
        for i=1:compsteps-jumps
            % --- reconstruction  and projection
            dyn = b.*exp(omega*time(i));
            X_rec(:,i) = Nu*dyn;

            % --- errors
            err_vec(i) = norm( X_comp(:,i+jumps)- X_rec(:,i) )/norm(X_comp(:,i+jumps));
              %index shift is necessary since the time domain version reconstructs
              % the 2. snapshot at t0
              %calculate the summed error for all calculated snaps 
            n_full = n_full + err_vec(i);
            count = count + 1;
            %calculate the summed error only for the reconstructed input snaps
            if count == jumps
                n = n + err_vec(i);
                count = 0;
            end 
        end
    end
end

%% calculate some evaluation parameters

if FULL == 1
    %sample variance
    samvar_full = std(err_vec);
    samvar = std(err_vec(1:jumps:end));
    %maximum
    maxis_full = max(err_vec);
    maxis = max(err_vec(1:jumps:end));
    %summed interpolation error
    n_int = n_full - n;
    %average errors on full interval
    nq_all = n_full/length(err_vec);
    nq = n/length(err_vec(1:jumps:end));
    %average errors on reconstruction interval from 0 to cycles
    n_rec_all = sum(err_vec(1:1:jumps*floor(cycles*50/jumps-1)));
    n_rec_rec = sum(err_vec(1:jumps:jumps*floor(cycles*50/jumps-1)));
    n_rec_int = n_rec_all - n_rec_rec;
    nq_rec_all = n_rec_all/(jumps*floor(cycles*50/jumps-1));
    nq_rec_rec = n_rec_rec/floor(cycles*50/jumps-1);
    nq_rec_int = n_rec_int/((jumps-1)*floor(cycles*50/jumps-1));


elseif FULL == 0
    %sample variance
    samvar_full = std(err_vec);
    %maximum
    maxis_full = max(err_vec);
    %average error
    nq_full = n/length(err_vec);
    %average errors on reconstruction interval from 0 to cycles
    n_rec = sum(err_vec(1:1:floor(cycles*50/jumps-1)));
    nq_rec = n_rec/(floor(cycles*50/jumps-1));

end

%count all unstable eigenvalues (with tolerance)
EWs=0;
for j = 1:length(diag(Lamb)) 
    if norm(Lamb(j,j)) > 1+1e-12
        EWs = EWs+1;
    end
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load input data
in_data = readmatrix('STAR_massFlowRate_schindler_5zyklen.txt');
[ ~, idx]= min( abs( in_data(:,1) - steps*dt ) );


%rel. error plots
figure()
hold on
yyaxis left
plot(in_data(1:idx,1),in_data(1:idx,2))
ylabel('MassFlowRate Input  $\dot{m}$', 'interpreter','latex', 'fontsize',16)
hold off
hold on
yyaxis right
if FULL == 0
    semilogy([dt+dt*jumps:dt*jumps:steps/spc],err_vec)
elseif FULL == 1
    semilogy([dt+dt*jumps:dt*jumps:steps/spc],err_vec(1:jumps:end))
end
semilogy([0 5],[1e-1 1e-1])
semilogy([0 5],[3e-1 3e-1], 'r')
hold off
set(gca,'yscale','log')
ylabel(sprintf('rel. error for $\\Delta t=%.2f$ , rank=%i',dt*jumps,rank), 'interpreter','latex', 'fontsize',16)
xlabel('time $t$ in $s$', 'interpreter','latex', 'fontsize',18)
title('DMD with exp. time projection', 'fontsize',24)
grid minor
h=gcf;
h.WindowState = 'maximized';
path = fullfile('D:\Eigene Dokumente\Uni\BA\',sprintf('%iaus%.1f_r%i_%i_%i_%i.fig',jumps,cycles,rank,ROM,FULL,FULL_Plot));
savefig(path)
close figure 1


if FULL_Plot == 1
    %rel. error plots with full interval additionally
    figure() 
    hold on
    yyaxis left
    plot(in_data(1:idx,1),in_data(1:idx,2))
    ylabel('MassFlowRate Input  $\dot{m}$', 'interpreter','latex', 'fontsize',16)
    hold off
    hold on
    yyaxis right
    semilogy([dt+dt*jumps:dt:steps/spc],err_vec)
    semilogy([dt+dt*jumps:dt*jumps:steps/spc],err_vec(1:jumps:end), 'Color','g','marker','*')
    semilogy([0 5],[1e-1 1e-1])
    semilogy([0 5],[3e-1 3e-1], 'r')
    hold off
    set(gca,'yscale','log')
    ylabel(sprintf('rel. error for $\\Delta t=%.2f$ , rank=%i, cycles=%.1f',dt,rank,cycles), 'interpreter','latex', 'fontsize',16)
    xlabel('time $t$ in $s$', 'interpreter','latex', 'fontsize',18)
    title('DMD with exp. time projection on full interval', 'fontsize',24)
    grid minor
    h=gcf;
    h.WindowState = 'maximized';
    path = fullfile('D:\Eigene Dokumente\Uni\BA\',sprintf('full_%iaus%.1f_r%i_%i_%i_%i.fig',jumps,cycles,rank,ROM,FULL,FULL_Plot));
    savefig(path)
    close figure 1
end



%eigenvalue plots
figure()
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'-') % plot unit circle
xlabel('Real($\Lambda(A)$)', 'interpreter','latex', 'fontsize',14)
ylabel('Imag($\Lambda(A)$)', 'interpreter','latex', 'fontsize',14)
title('eigenvalues for DMD set','fontsize',17)
hold on, grid on
scatter(real(diag(Lamb)), imag(diag(Lamb)),'o')
axis([-1.1 1.1 -1.1 1.1]);
annotation('textbox', [0.78,0.79, 0.15, 0.15], 'String', sprintf("rank = %i\njumps = %i\ncycles = %.1f",rank,jumps,cycles),'FitBoxToText','on','FontSize',10,'BackgroundColor','w')
path = fullfile('D:\Eigene Dokumente\Uni\BA\',sprintf('EW_%iaus%.1f_r%i_%i_%i_%i.fig',jumps,cycles,rank,ROM,FULL,FULL_Plot));
savefig(path)
close figure 1

 
end