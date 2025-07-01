clear; close all

%% FIGURE PARAMETERS

set( 0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0 );

%% MFS

% define the grid 

N = 512;           % must be even 
L = 400;
delta = 0.3;
dx = L / N;
x = - L / 2 : dx : L / 2 - dx;
x = x';  

%% PHYSICAL PARAMETERS

% wavenumbers

Nf = 41;
f = linspace( 3.1, 5.1, Nf );   % 3.1 GHz - 5.1 GHz
c0 = 29.9792458;                % speed of light in centimeters GHz
k0 = 2.0 * pi * f / c0;         % wavenumber for z > zb

% synthetic aperture

Na = 35;
a = 102;                           % aperture
xa = linspace( -a/2, a/2, Na );    % path locations
za = 75;                           % absolute height of synthetic aperture

% reflectivity

rho = 1;

% location of target - "point-like"

y0_1 = 25;
y0_2 = 10;
z0_values = linspace(5, 50, 5);

% allocate space for measurements

Dss = zeros( Nf, Na );
Dsh = zeros( Nf, Na );

Rss = zeros( Nf, Na );
Rsh = zeros( Nf, Na );

S = zeros( Nf, Na );

Sss = zeros( Nf, Na );
Ssh = zeros( Nf, Na );

Sigss = zeros( Nf, Na );
Sigsh = zeros( Nf, Na );

Uss = zeros( Nf, Na );
Ush = zeros( Nf, Na );

Vss = zeros( Nf, Na );
Vsh = zeros( Nf, Na );

Dss_flat = zeros( Nf, Na );
Dsh_flat = zeros( Nf, Na );

Rss_flat = zeros( Nf, Na );
Rsh_flat = zeros( Nf, Na );

Sss_flat = zeros( Nf, Na );
Ssh_flat = zeros( Nf, Na );

exact_ss = zeros( Nf, Na );
exact_sh = zeros( Nf, Na );

% pre-computed distances 

[ ii, jj ] = ndgrid( (1:N) );              % ii: surface points, jj: source points
[ ia, ja ] = ndgrid( (1:N), (1:Na) );

%% ROUGH SURFACE

rng( 'default' );

hRMS_values = [ 1e-3 1e-2 1e-1 5e-1 1 ]; %linspace(0.1, 0.5, 5);        % 0 for flat surface 
delta_Rss = zeros(size(hRMS_values));
delta_Rsh = zeros(size(hRMS_values));
delta_S = zeros(size(hRMS_values));

sigs_ss = zeros(2,length(hRMS_values));
sigs_sh = zeros(2,length(hRMS_values));

for i = 1 : length(hRMS_values)

    hRMS = hRMS_values(i);
    z0 = z0_values(i);
    lcor = 8.0;

    % correlation function

    C = exp( - 0.5 * x.^2 / lcor^2 );

    % zero mean circularly symmetric complex Gaussian random variables

    zeta = ( rand(1,N/2-1) + 1j * randn(1,N/2-1) ) / sqrt( 2 );
    zeta = [ randn(1) zeta randn(1) fliplr(conj(zeta) ) ]';

    % random function

    z = real( ifft( sqrt( fft( C ) ) .* zeta ) );
    mu = mean( z )';
    sig = std( z )';
    z = hRMS * ( z - mu ) / sig;

    %% DEFINING THE PROBLEM

    % derivatives

    xk = 2.0 * pi / L * fftshift( -N/2 : N/2-1 );
    xk = xk';
    dz = real( ifft( 1j * xk .* fft( z ) ) );

    % components of the unit normal

    nu_x = -dz ./ sqrt( 1 + dz.^2 );
    nu_z = 1 ./ sqrt( 1 + dz.^2 );

    %% COMPUTE MEASUREMENTS

    % pre-compute distances for rough surface

    R = sqrt( ( x(ii) - x(jj) ).^2 + ( z(ii) - z(jj) + delta ).^2 );

    R_flat = sqrt( ( x(ii) - x(jj) ).^2 + delta^2 );

    A = sqrt(( x(ia) - xa(ja) ).^2 + ( z(ia) - za ).^2);

    A_flat = sqrt(( x(ia) - xa(ja) ).^2 + za^2);

    % loop over frequencies and locations

    for m = 1 : Nf

        T = sqrt((xa(ja) - y0_1).^2 + (za - z0).^2);

        % first block - reflected field

        H = 1j / 4 * besselh( 0, 1, k0(m) * R );

        H_flat = 1j / 4 * besselh( 0, 1, k0(m) * R_flat );

        H_prime = ( -1j * k0(m) / 4  * nu_x(ii) .* (  x(ii) - x(jj)  ) ./ R ...
            .* besselh( 1, 1, k0(m) * R ) ) + ( -1j * k0(m) / 4 * nu_z(ii) ...
            .* (  z(ii) - z(jj) + delta ) ./ R .* besselh( 1, 1, k0(m) * R ) );

        H_prime_flat = ( -1j * k0(m) / 4 * ...
            ( delta ) ./ R_flat .* besselh( 1, 1, k0(m) * R_flat ) );

        % second block - scattered field

        J = sqrt( ( x - y0_1 ).^2 + ( z - y0_2 ).^2 );

        J_flat = sqrt( ( x - y0_1 ).^2 + y0_2^2 );

        b = 1j / 4 * besselh( 0, 1, k0(m) * J );

        b_flat = 1j / 4 * besselh( 0, 1, k0(m) * J_flat );

        db = -1j * k0(m) / 4 * nu_x .* ( x - y0_1 ) ./ J .* besselh( 1, 1, k0(m) * J ) ...
            -1j * k0(m) / 4 * nu_z .* ( z - y0_2 ) ./ J .* besselh( 1, 1, k0(m) * J );

        db_flat = -1j * k0(m) / 4 .* ( -y0_2 ) ./ J_flat .* besselh( 1, 1, k0(m) * J_flat );

        % third block - self-consistent argument

        P = sqrt( ( y0_1 - x ).^2 + ( y0_2 - z + delta ).^2 );

        P_flat = sqrt( ( y0_1 - x ).^2 + ( y0_2 + delta ).^2 );

        c = 1j / 4 * besselh( 0, 1, k0(m) * P );    % sound-soft and sound-hard same;
                                                    % (not evaluated on boundary)
        c_flat = 1j / 4 * besselh( 0, 1, k0(m) * P_flat );

        c = c.';    % want row vector
        c_flat = c_flat.';

        % right-hand side

        M = sqrt( ( xa - y0_1 ).^2 + ( za - y0_2 ).^2 );
        u_i_y0 = 1j / 4 * besselh( 0, 1, k0(m) * M );   % sound-soft = sound-hard

        % incident field

        u_i = 1j / 4 * besselh( 0, 1, k0(m) * A );

        u_i_flat = 1j / 4 * besselh( 0, 1, k0(m) * A_flat );

        du_i = -1j * k0(m) / 4 * nu_x .* ( x(ia) - xa(ja) ) ./ A .* besselh( 1, 1, k0(m) * A ) ...
            -1j * k0(m) / 4 * nu_z .* ( z(ia) - za ) ./ A .* besselh( 1, 1, k0(m) * A );

        du_i_flat = -1j * k0(m) / 4 .* ( -za ) ./ A_flat .* besselh( 1, 1, k0(m) * A_flat );

        % sound-soft rough
        ss_matrix = [ H  rho * b  ;  c  -1 ];
        ss_rhs = [ -u_i ; -u_i_y0 ];

        w_ss = ss_matrix \ ss_rhs;
        u_e_ss = w_ss(end,:);      % exciting field

        % sound-soft flat
        ss_matrix_flat = [ H_flat  rho * b_flat  ;  c_flat  -1 ];
        ss_rhs_flat = [ -u_i_flat ; -u_i_y0 ];

        w_ss_flat = ss_matrix_flat \ ss_rhs_flat;

        % sound-hard rough
        sh_matrix = [ H_prime  rho * db  ;  c  -1 ];
        sh_rhs = [ -du_i ; -u_i_y0 ];

        w_sh = sh_matrix \ sh_rhs;
        u_e_sh = w_sh(end,:);

        % sound-hard flat
        sh_matrix_flat = [ H_prime_flat  rho * db_flat  ;  c_flat  -1 ];
        sh_rhs_flat = [ -du_i_flat ; -u_i_y0 ];

        w_sh_flat = sh_matrix_flat \ sh_rhs_flat;

        % evaluate solution for measurements

        % sound-soft rough
        Sss(m,:) = diag( rho * 1j / 4 * besselh(0, 1, k0(m) * T) .* u_e_ss(1,:) );     % scattered
        Rss(m,:) = diag( 1j / 4 * besselh( 0, 1, k0(m) * A ).' * w_ss(1:end-1,:) );    % reflected
        Dss(m,:) = Rss(m,:) + Sss(m,:);                                                % full measurement

        % sound-soft flat
        Sss_flat(m,:) = diag( rho * 1j / 4 * besselh(0, 1, k0(m) * T ) .* u_e_ss(1,:) );                % scattered
        Rss_flat(m,:) = diag( 1j / 4 * besselh( 0, 1, k0(m) * A_flat ).' * w_ss_flat(1:end-1,:) );      % reflected
        Dss_flat(m,:) = Rss_flat(m,:) + Sss_flat(m,:);                                                  % full measurement

        % sound-hard rough
        Ssh(m,:) = diag( rho * 1j / 4 * besselh(0, 1, k0(m) * T ) .* u_e_sh(1,:) );
        Rsh(m,:) = diag( 1j / 4 * besselh( 0, 1, k0(m) * A ).' * w_sh(1:end-1,:) );
        Dsh(m,:) = Rsh(m,:) + Ssh(m,:);

        % sound-hard flat 
        Ssh_flat(m,:) = diag( rho * 1j / 4 * besselh(0, 1, k0(m) * T) .* u_e_sh(1,:) );
        Rsh_flat(m,:) = diag( 1j / 4 * besselh( 0, 1, k0(m) * A_flat ).' * w_sh_flat(1:end-1,:) );
        Dsh_flat(m,:) = Rsh_flat(m,:) + Ssh_flat(m,:);

        % analytic solutions (for flat surface)
        exact_ss(m,:) = - 1j / 4 ...
            * besselh( 0, 1, k0(m) * 2 * za ) ;

        exact_sh(m,:) = 1j / 4 ...
            * besselh( 0, 1, k0(m) * 2 * za );

        % matrix of SAR measurements in the whole space due to a point target (y0)
        S(m,:) = diag( rho * (1j / 4 * besselh(0, 1, k0(m) * T )).^2 );

     end

    %% ERROR ANALYSIS

    RE_ss = abs( Rss -  Rss_flat );
    RE_sh = abs( Rsh -  Rsh_flat );

    delta_Rss(i) = norm( RE_ss, "fro" );
    delta_Rsh(i) = norm( RE_sh, "fro" );

    delta_S(i) = norm( S, "fro" );

    %% SINGULAR VALUES

    [ ~, Sigss, ~ ] = svd( Dss );
    [ ~, Sigsh, ~ ] = svd( Dsh );
 
    sigs_ss(:,i) = diag(Sigss(1:2,1:2));
    sigs_sh(:,i) = diag(Sigsh(1:2,1:2));

end

%% PLOT SINGULAR VALUES VS. HRMS

figure(1);
loglog( hRMS_values, sigs_ss, '-o', 'LineWidth', 2);
xlabel('$$\log(h_{RMS})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
ylabel('$$\log(\sigma^{ss})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
title('$$\sigma^{ss}$$ vs Surface Roughness', 'Interpreter', 'LaTeX', 'fontsize', 24); 
hold on
% calculate best fit line and plot 
p2 = polyfit( log(hRMS_values(3:5)), log(sigs_ss(2,3:5)), 1);
slope = p2(1);
loglog(hRMS_values(3:5), exp(polyval(p2, log(hRMS_values(3:5)))), '--r', 'LineWidth', 2);
slope_label = sprintf('Slope: %.2f', slope);
annotation('textbox', [0.35, 0.3, 0.2, 0.1], 'String', slope_label, 'BackgroundColor', 'none', 'edgecolor', 'none', 'fontsize', 16, 'fontname', 'times');
legend('$$\sigma_{1}^{ss}$$', '$$\sigma_{2}^{ss}$$', 'Best Fit Line','Interpreter', 'LaTeX', 'Location', 'best')
hold off

figure(2);
loglog( hRMS_values, sigs_sh, '-o', 'LineWidth', 2);
xlabel('$$\log(h_{RMS})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
ylabel('$$\log(\sigma^{sh})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
title('$$\sigma^{sh}$$ vs Surface Roughness', 'Interpreter', 'LaTeX', 'fontsize', 24);
hold on
% calculate best fit line and plot 
p2 = polyfit( log(hRMS_values(3:5)), log(sigs_sh(2,3:5)), 1);
slope = p2(1);
loglog(hRMS_values(3:5), exp(polyval(p2, log(hRMS_values(3:5)))), '--r', 'LineWidth', 2);
slope_label = sprintf('Slope: %.2f', slope);
annotation('textbox', [0.35, 0.3, 0.2, 0.1], 'String', slope_label, 'BackgroundColor', 'none', 'edgecolor', 'none', 'fontsize', 16, 'fontname', 'times');
legend('$$\sigma_{1}^{sh}$$', '$$\sigma_{2}^{sh}$$', 'Best Fit Line', 'Interpreter', 'LaTeX', 'Location', 'best')
hold off

%% SUBTRACT FLAT SURFACE MEASUREMENTS

Dbar_ss = Dss - Rss_flat;
Dbar_sh = Dsh - Rsh_flat;

D0_ss = Dss - exact_ss;
D0_sh = Dsh - exact_sh;

%% KIRCHHOFF MIGRATION

Nxgrid = 101;
Nzgrid = 101;

x_grid = linspace(-50,50,Nxgrid);
z_grid = linspace(1,50,Nzgrid);

% mesh grid
[ Xmesh, Zmesh ] = meshgrid( x_grid, z_grid );

% make space for image

KM_Rss = 0 * Xmesh;
KM_Rsh = 0 * Xmesh;

KM_Sss = 0 * Xmesh;
KM_Ssh = 0 * Xmesh;

KM_Dss = 0 * Xmesh;
KM_Dsh = 0 * Xmesh;

KM_Dbar_ss = 0 * Xmesh;
KM_Dbar_sh = 0 * Xmesh;

KM_D0_ss = 0 * Xmesh;
KM_D0_sh = 0 * Xmesh;

for m = 1 : Nf

    for n = 1 : Na

        Rn = sqrt( ( Xmesh - xa(n) ).^2 + ( Zmesh - za ).^2 );

        KM_Rss = KM_Rss + Rss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_Rsh = KM_Rsh + Rsh(m,n) * exp( -1j * 2 * k0(m) * Rn );

        KM_Sss = KM_Sss + Sss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_Ssh = KM_Ssh + Ssh(m,n) * exp( -1j * 2 * k0(m) * Rn );

        KM_Dss = KM_Dss + Dss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_Dsh = KM_Dsh + Dsh(m,n) * exp( -1j * 2 * k0(m) * Rn );

        KM_Dbar_ss = KM_Dbar_ss + Dbar_ss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_Dbar_sh = KM_Dbar_sh + Dbar_sh(m,n) * exp( -1j * 2 * k0(m) * Rn );

        KM_D0_ss = KM_D0_ss + D0_ss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_D0_sh = KM_D0_sh + D0_sh(m,n) * exp( -1j * 2 * k0(m) * Rn );

    end

end


%% PLOT THE IMAGES

figure(3)
pcolor(x_grid, z_grid, abs(KM_Dss));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('KM on D (Sound-Soft)', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

figure(4)
pcolor(x_grid, z_grid, abs(KM_Dsh));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('KM on D (Sound-Hard)', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

% with R_flat subtracted 

figure(5)
pcolor(x_grid, z_grid, abs(KM_Dbar_ss));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('KM on $$\bar{D}$$ (Sound-Soft)', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

figure(6)
pcolor(x_grid, z_grid, abs(KM_Dbar_sh));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('KM on $$\bar{D}$$ (Sound-Hard)', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

% with analytic solution subtracted 

figure(7)
pcolor(x_grid, z_grid, abs(KM_D0_ss));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('KM on $$D_0$$ (Sound-Soft)', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

figure(8)
pcolor(x_grid, z_grid, abs(KM_D0_sh));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('KM on $$D_0$$ (Sound-Hard)', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

%% PCA/SVD TO REMOVE REFLECTIONS

[ Uss, Sigss, Vss ] = svd( Dss );
[ Ush, Sigsh, Vsh ] = svd( Dsh );

% plot the singular values

figure(9)
semilogy( (1:min(Nf,Na)), diag(Sigss), 'o', (1:min(Nf,Na)), diag(Sigsh), 'x' )
grid on
xlabel( '$j$', 'Interpreter', 'LaTeX', 'fontsize', 24 );
ylabel( '$\sigma_{j}$', 'Interpreter', 'LaTeX', 'fontsize', 24 );
legend( {'sound-soft', 'sound-hard' }, 'Interpreter', 'LaTeX', 'fontsize', 24 );
title( 'Singular Values of $D$', 'Interpreter', 'LaTeX', 'fontsize', 24 )

% set the truncation level

ntrunc = 3;

Dtilde_ss = Dss;
Dtilde_sh = Dsh;

for j = 1 : ntrunc

    Dtilde_ss = Dtilde_ss - Sigss(j,j) * Uss(:,j) * Vss(:,j)';
    Dtilde_sh = Dtilde_sh - Sigsh(j,j) * Ush(:,j) * Vsh(:,j)';

end

%% KM ON PCA PROCESSED DATA

% make space for image

KM_DPCAss = 0 * Xmesh;
KM_DPCAsh = 0 * Xmesh;

for m = 1 : Nf

    for n = 1 : Na

        Rn = sqrt( ( Xmesh - xa(n) ).^2 + ( Zmesh - za ).^2 );

        KM_DPCAss = KM_DPCAss + Dtilde_ss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_DPCAsh = KM_DPCAsh + Dtilde_sh(m,n) * exp( -1j * 2 * k0(m) * Rn );

    end

end

figure(10)

pcolor(x_grid, z_grid, abs(KM_DPCAss));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15)
hold off
title('PCA sound soft', 'Interpreter', 'LaTeX', 'fontsize', 24);
subtitle('$\sigma_j$ truncated: 2', 'Interpreter', 'LaTeX', 'fontsize', 18)
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

figure(11)

pcolor(x_grid, z_grid, abs(KM_DPCAsh));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('PCA sound hard', 'Interpreter', 'LaTeX', 'fontsize', 24);
subtitle('$\sigma_j$ truncated: 2', 'Interpreter', 'LaTeX', 'fontsize', 18)
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

%% ERROR PLOT 

figure(12)

loglog( hRMS_values, delta_Rss, '-o', 'LineWidth', 2);
xlabel('$$\log(h_{RMS})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
ylabel('$$\log(\delta_{Rss})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
title('Error vs Surface Roughness', 'Interpreter', 'LaTeX', 'fontsize', 24);
hold on
% calculate best fit line and plot 
p = polyfit( log(hRMS_values), log(delta_Rss), 1);
slope = p(1);
loglog(hRMS_values, exp(polyval(p, log(hRMS_values))), '--r', 'LineWidth', 2);
slope_label = sprintf('Slope: %.2f', slope);
annotation('textbox', [0.35, 0.6, 0.2, 0.1], 'String', slope_label, 'BackgroundColor', 'none', 'edgecolor', 'none', 'fontsize', 16, 'fontname', 'times');
hold off

figure(13)

loglog( hRMS_values, delta_Rsh, '-o', 'LineWidth', 2);
xlabel('$$\log(h_{RMS})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
ylabel('$$\log(\delta_{Rsh})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
title('Error vs Surface Roughness', 'Interpreter', 'LaTeX', 'fontsize', 24); 
hold on
% calculate best fit line and plot 
p2 = polyfit( log(hRMS_values), log(delta_Rsh), 1);
slope = p2(1);
loglog(hRMS_values, exp(polyval(p2, log(hRMS_values))), '--r', 'LineWidth', 2);
slope_label = sprintf('Slope: %.2f', slope);
annotation('textbox', [0.35, 0.6, 0.2, 0.1], 'String', slope_label, 'BackgroundColor', 'none', 'edgecolor', 'none', 'fontsize', 16, 'fontname', 'times');
hold off

figure(14)
loglog( abs(za - z0_values), delta_S, '-o', 'LineWidth', 2);
xlabel('$$\log(za - z0)$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
ylabel('$$\log(\delta_{S})$$', 'Interpreter', 'LaTeX', 'fontsize', 24);
title('Error vs Target Location', 'Interpreter', 'LaTeX', 'fontsize', 24);
hold on
% calculate best fit line and plot 
p = polyfit( log(abs(za - z0_values)), log(delta_S), 1);
slope = p(1);
loglog( abs(za - z0_values), exp(polyval(p, log(abs(za - z0_values)))), '--r', 'LineWidth', 2);
slope_label = sprintf('Slope: %.2f', slope);
annotation('textbox', [0.45, 0.6, 0.2, 0.1], 'String', slope_label, 'BackgroundColor', 'none', 'edgecolor', 'none', 'fontsize', 16, 'fontname', 'times');
hold off

