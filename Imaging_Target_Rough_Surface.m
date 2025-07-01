clear; close all

% Imaging Final Project

%% FIGURE PARAMETERS

set( 0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0 );

%% MFS

% define the grid 
N = 512;           % must be even 
L = 400;
dx = L / N;
delta = 0.3;
x = - L / 2 : dx : L / 2 - dx;
x = x';  

%% PHYSICAL PARAMETERS

% wavenumbers
Nf   = 41;
f    = linspace( 3.1, 5.1, Nf );  % 3.1 GHz - 5.1 GHz
c0   = 29.9792458;                % speed of light in centimeters GHz
k0   = 2.0 * pi * f / c0;         % wavenumber for z > zb

% synthetic aperture
Na = 35;
a  = 102;                          % aperture
xa = linspace( -a/2, a/2, Na );    % path locations
za = 75;                           % absolute height of synthetic aperture

%% ROUGH SURFACE

rng( 'default' );

hRMS = 0.4;        % 0 for flat surface  

lcor = 6.0;

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

% reflectivity 
rho = 1; 

% location of target - "point-like"

y0_1 = 25;
y0_2 = 20;

%% COMPUTE MEASUREMENTS

% allocate space for measurements

Dss = zeros( Nf, Na );
Dsh = zeros( Nf, Na );

Rss = zeros( Nf, Na );
Rsh = zeros( Nf, Na );

Sss = zeros( Nf, Na );
Ssh = zeros( Nf, Na );

% pre-compute distances for surface

[ ii, jj ] = ndgrid( (1:N) );    % ii: surface points, jj: source points

R = sqrt( ( x(ii) - x(jj) ).^2 + ( z(ii) - z(jj) + delta ).^2 );

% pre-computed distances to measurements

[ ia, ja ] = ndgrid( (1:N), (1:Na) );

A = ( x(ia) - xa(ja) ).^2 + ( z(ia) - za ).^2;

% loop over frequencies and locations

for m = 1 : Nf

    % first block - reflected field
    
    H = 1j / 4 * besselh( 0, 1, k0(m) * R );

    H_prime = ( -1j * k0(m) / 4  * nu_x(ii) .* (  x(ii) - x(jj)  ) ./ R ...
        .* besselh( 1, 1, k0(m) * R ) ) + ( -1j * k0(m) / 4 * nu_z(ii) ...
        .* (  z(ii) - z(jj) + delta ) ./ R .* besselh( 1, 1, k0(m) * R ) );

    % second block - scattered field

    J = sqrt( ( x - y0_1 ).^2 + ( z - y0_2 ).^2 );

    b = 1j / 4 * besselh( 0, 1, k0(m) * J );  

    db = -1j * k0(m) / 4 * nu_x .* ( x - y0_1 ) ./ J .* besselh( 1, 1, k0(m) * J ) ...
        -1j * k0(m) / 4 * nu_z .* ( z - y0_2 ) ./ J .* besselh( 1, 1, k0(m) * J );

    % third block - self-consistent argument

    P = sqrt( ( y0_1 - x ).^2 + ( y0_2 - z + delta ).^2 );

    c = 1j / 4 * besselh( 0, 1, k0(m) * P );    % sound-soft and sound-hard same; 
                                                % (not evaluated on boundary) 
    c = c.';    % want row vector 

    % right-hand side 

    M = sqrt( ( xa - y0_1 ).^2 + ( za - y0_2 ).^2 );
    u_i_y0 = 1j / 4 * besselh( 0, 1, k0(m) * M );   % sound-soft = sound-hard

    % incident field

    u_i = 1j / 4 * besselh( 0, 1, k0(m) * sqrt( A ) );

    du_i = -1j * k0(m) / 4 * nu_x .* ( x(ia) - xa(ja) ) ./ sqrt(A) .* besselh( 1, 1, k0(m) * sqrt(A) ) ...
        -1j * k0(m) / 4 * nu_z .* ( z(ia) - za ) ./ sqrt(A) .* besselh( 1, 1, k0(m) * sqrt(A) );

    % sound-soft
    ss_matrix = [ H  rho * b  ;  c  -1 ];
    ss_rhs = [ -u_i ; -u_i_y0 ];

    w_ss = ss_matrix \ ss_rhs;
    u_e_ss = w_ss(end,:);      % exciting field 

    % sound-hard
    sh_matrix = [ H_prime  rho * db  ;  c  -1 ];
    sh_rhs = [ -du_i ; -u_i_y0 ];

    w_sh = sh_matrix \ sh_rhs;
    u_e_sh = w_sh(end,:); 

    % evaluate solution for measurements

    for n = 1 : Na

        K = sqrt((x - xa(n)).^2 + (z - za).^2);
        T = sqrt((xa(n) - y0_1).^2 + (za - y0_2).^2);

        % sound-soft
        Rss(m,n) = 1j / 4 * besselh(0, 1, k0(m) * K).' * w_ss(1:end-1,n);   % reflected
        Sss(m,n) = rho * 1j / 4 * besselh(0, 1, k0(m) * T) * u_e_ss(n);     % scattered
        Dss(m,n) = Rss(m,n) + Sss(m,n);                                     % full measurement

        % sound-hard
        Rsh(m,n) = 1j / 4 * besselh(0, 1, k0(m) * K).' * w_sh(1:end-1,n);
        Ssh(m,n) = rho * 1j / 4 * besselh(0, 1, k0(m) * T) * u_e_sh(n);
        Dsh(m,n) = Rsh(m,n) + Ssh(m,n);

    end

end

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

for m = 1 : Nf

    for n = 1 : Na

        Rn = sqrt( ( Xmesh - xa(n) ).^2 + ( Zmesh - za ).^2 );
        
        KM_Rss = KM_Rss + Rss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_Rsh = KM_Rsh + Rsh(m,n) * exp( -1j * 2 * k0(m) * Rn );

        KM_Sss = KM_Sss + Sss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_Ssh = KM_Ssh + Ssh(m,n) * exp( -1j * 2 * k0(m) * Rn );

        KM_Dss = KM_Dss + Dss(m,n) * exp( -1j * 2 * k0(m) * Rn );
        KM_Dsh = KM_Dsh + Dsh(m,n) * exp( -1j * 2 * k0(m) * Rn );

    end

end

%% PLOT THE IMAGES
% 
% figure(1)
% pcolor(x_grid, z_grid, abs(KM_Rss));
% shading flat;
% hold on
% plot( y0_1, y0_2, 'r*', 'markersize', 15);
% hold off
% title('Reflection', 'Interpreter', 'LaTeX', 'fontsize', 15);
% xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 15); 
% ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 15);
% colorbar;

figure(2)
pcolor(x_grid, z_grid, abs(KM_Sss));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('Target Scattering: sound-soft', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

figure(3)
pcolor(x_grid, z_grid, abs(KM_Dss));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('Full Measurement: sound-soft', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

% figure(4)
% pcolor(x_grid, z_grid, abs(KM_Rsh));
% shading flat;
% hold on
% plot( y0_1, y0_2, 'r*', 'markersize', 15);
% hold off
% title('Reflection', 'Interpreter', 'LaTeX', 'fontsize', 15);
% xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 15); 
% ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 15);
% colorbar;

figure(5)
pcolor(x_grid, z_grid, abs(KM_Ssh));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('Target Scattering: sound-hard', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

figure(6)
pcolor(x_grid, z_grid, abs(KM_Dsh));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('Full Measurement: sound-hard', 'Interpreter', 'LaTeX', 'fontsize', 24);
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

%% PCA/SVD TO REMOVE REFLECTIONS

[ Uss, Sigss, Vss ] = svd( Dss );
[ Ush, Sigsh, Vsh ] = svd( Dsh );

% plot the singular values

figure(7)
semilogy( (1:min(Nf,Na)), diag(Sigss), 'o', (1:min(Nf,Na)), diag(Sigsh), 'x' )
grid on
xlabel( '$j$', 'Interpreter', 'LaTeX', 'fontsize', 24 );
ylabel( '$\sigma_{j}$', 'Interpreter', 'LaTeX', 'fontsize', 24 );
legend( {'sound-soft', 'sound-hard' }, 'Interpreter', 'LaTeX', 'fontsize', 24 );
title( 'Singular Values of $D$', 'Interpreter', 'LaTeX', 'fontsize', 24 )

% set the truncation level

ntrunc = 1;

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

figure(8)

pcolor(x_grid, z_grid, abs(KM_DPCAss));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15)
hold off
title('Sound-Soft: PCA', 'Interpreter', 'LaTeX', 'fontsize', 24);
subtitle('$\sigma_j$ truncated: 1', 'Interpreter', 'LaTeX', 'fontsize', 18)
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;

figure(9)

pcolor(x_grid, z_grid, abs(KM_DPCAsh));
shading flat;
hold on
plot( y0_1, y0_2, 'r*', 'markersize', 15);
hold off
title('Sound-Hard: PCA', 'Interpreter', 'LaTeX', 'fontsize', 24);
subtitle('$\sigma_j$ truncated: 1', 'Interpreter', 'LaTeX', 'fontsize', 18)
xlabel('$x$', 'Interpreter', 'LaTeX', 'fontsize', 24); 
ylabel('$z$', 'Interpreter', 'LaTeX', 'fontsize', 24);
colorbar;
