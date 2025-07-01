clear; close all

% Method of Fundamental Solutions for Rough Surface Scattering 

%% FIGURE PARAMETERS

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0);

%% MFS

% define the grid 
N = 300;           % must be even 
L = 100;
dx = L / N;        % delta - small perturbation
x = - L / 2 : dx : L / 2 - dx;
x = x';            % want column vector 

%% ROUGH SURFACE

rng( 'shuffle' );
hRMS = 0.9;        % 0 for flat surface  

lcor = 8;

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

% derivatives
xk = 2.0 * pi / L * fftshift( -N/2 : N/2-1 );
xk = xk';
dz = real( ifft( 1j * xk .* fft( z ) ) );

% components of the unit normal
nu_x = -dz ./ sqrt( 1 + dz.^2 );
nu_z = 1 ./ sqrt( 1 + dz.^2 );

% wave constant 
k = 1;

% initial source point
x0 = 20; 
z0 = 25;

% compute the incident field & derivative (deriv for SH)
A = ( x - x0 ).^2 + ( z - z0 ).^2;
u_i = 1j / 4 * besselh( 0, 1, k * sqrt( A ) );
du_i = -1j * k / 4 * nu_x .* ( x - x0 ) ./ sqrt(A) .* besselh( 1, 1, k * sqrt(A) ) ...
    -1j * k / 4 * nu_z .* ( z - z0 ) ./ sqrt(A) .* besselh( 1, 1, k * sqrt(A) );

% compute the MFS system matrix & derivative (deriv for SH)
[ ii, jj ] = ndgrid( (1:N) );

R = sqrt( ( x(ii) - x(jj) ).^2 + ( z(ii) - z(jj) + dx ).^2 );

H = 1j / 4 * besselh( 0, 1, k * R );

H_prime = ( -1j * k / 4  * nu_x(ii) .* (  x(ii) - x(jj)  ) ./ R ...
    .* besselh( 1, 1, k * R ) ) + ( -1j * k / 4 * nu_z(ii) ...
    .* (  z(ii) - z(jj) + dx ) ./ R .* besselh( 1, 1, k * R ) );

% flat surface matrices
R0 = sqrt( ( x(ii) - x(jj) ).^2 + dx^2 );

H0 = 1j / 4 * besselh( 0, 1, k * R0 );

H0_prime = -1j * k / 4 * dx ./ R0 .* besselh( 1, 1, k * R0 );

test = 1j / 4 * sqrt( 2 ./ ( pi * k * R0 ) ) * exp( 1j * k * R0 + (1j*pi/4));   % behavior of H as z -> infinity

%% SOUND-SOFT BC

% weights are obtained from linear system -u_i = H*w
w_SS = H \ -u_i;

%% SOUND-HARD BC

% weights are obtained from linear system -du_prime = H_prime*w
w_SH = H_prime \ -du_i;

%% CREATE GRID POINTS FOR PLOT

Nxgrid = 101;
Nzgrid = 201;

x_grid = linspace(-50,50,Nxgrid);
z_grid = linspace(1,50,Nzgrid);

% mesh grid
[ Xmesh, Zmesh ] = meshgrid( x_grid, z_grid );

% useful index arrays
[ indx, jndx ] = ndgrid( (1:Nxgrid*Nzgrid), (1:N) );

% compute the reflected field over the mesh and normal derivative
R = sqrt( ( Xmesh(indx) - x(jndx) ).^2 + ( Zmesh(indx) - z(jndx) + dx ).^2 );

u_r_SS = 1j / 4 * besselh( 0, 1, k * R ) * w_SS;    % sound-soft
u_r_SH = 1j / 4 * besselh( 0, 1, k * R ) * w_SH;    % sound-hard

% reshape 
u_r_SS = reshape( u_r_SS, Nzgrid, Nxgrid );
u_r_SH = reshape( u_r_SH, Nzgrid, Nxgrid );

%% COMPUTE THE EXACT SOLUTIONS

% exact for flat surface 
u_r_exact_SS = - 1j / 4 ...
    * besselh( 0, 1, k * sqrt( ( Xmesh - x0 ).^2 + ( Zmesh + z0 ).^2 ) ); 

u_r_exact_SH = 1j / 4 ...
    * besselh( 0, 1, k * sqrt( ( Xmesh - x0 ).^2 + ( Zmesh + z0 ).^2 ) ); 

%% FIGURES 

% sound-soft bc
figure(1)
% subplot( 1, 3, 1 );
pcolor( x_grid, z_grid, real( u_r_exact_SS ) );
shading flat;
colorbar;      % gives height of surface
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( '$u_{flat}$: sound-soft', 'Interpreter', 'LaTeX' );
% sgtitle('Sound-soft', 'Interpreter', 'LaTeX', 'fontsize', 18)

figure(2)
% subplot( 1, 3, 2 );
pcolor( x_grid, z_grid, real( u_r_SS ) );
shading flat;
colorbar;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( '$u_r$: sound-soft', 'Interpreter', 'LaTeX' );

% subplot( 1, 3, 3 );
% pcolor( x_grid, z_grid, log10( abs( u_r_SS - u_r_exact_SS ) ) );
% shading flat;
% colorbar;
% hold on;
% plot( [-L/2 L/2 ], [ 1 1 ], 'o' )
% hold off;
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% title( 'Absolute Difference', 'Interpreter', 'LaTeX' );

% figure(2)
% 
% subplot(1,2,1)
% pcolor( x_grid, z_grid, real(u_r_SS) )
% shading flat;
% colorbar;
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% title( '2D (Sound-Soft)', 'Interpreter', 'LaTeX' );
% 
% subplot(1,2,2)
% surf( x_grid, z_grid, real(u_r_SS) )
% shading flat;
% colorbar;
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% title( '3D (Sound-Soft)', 'Interpreter', 'LaTeX' );

%sound-hard bc
figure(3)
% subplot( 1, 3, 1 );
pcolor( x_grid, z_grid, real( u_r_exact_SH ) );
shading flat;
colorbar;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( '$u_{flat}$: sound-hard', 'Interpreter', 'LaTeX' );
% sgtitle('Sound-hard', 'Interpreter', 'LaTeX', 'fontsize', 18)
% 
% subplot( 1, 3, 2 );
figure(4)
pcolor( x_grid, z_grid, real( u_r_SH ) );
shading flat;
colorbar;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( '$u_r$: sound-hard', 'Interpreter', 'LaTeX' );
% 
% subplot( 1, 3, 3 );
% pcolor( x_grid, z_grid, log10( abs( u_r_SH - u_r_exact_SH ) ) );
% shading flat;
% colorbar;
% hold on;
% plot( [ -L/2 L/2 ], [ 1 1 ], 'o' )
% hold off;
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% title( 'Absolute Difference', 'Interpreter', 'LaTeX' );
% 
% figure(4)
% 
% subplot(1,2,1)
% pcolor( x_grid, z_grid, real(u_r_SH) )
% shading flat;
% colorbar;
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% title( '2D (Sound-Hard)', 'Interpreter', 'LaTeX' );

% subplot(1,2,2)
% surf( x_grid, z_grid, real(u_r_SH) )
% shading flat;
% colorbar;
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% zlabel( '$u(x,z)$', 'Interpreter', 'LaTeX' );
% title( '3D (Sound-Hard)', 'Interpreter', 'LaTeX' );

% % plot of H and H_prime 
% figure(5)
% 
% subplot(1,2,1)
% imagesc( abs(H) );
% colorbar;
% title( '$H$', 'Interpreter', 'LaTeX' );
% 
% subplot(1,2,2)
% imagesc( abs(H_prime) );
% colorbar;
% title( '$H Prime$', 'Interpreter', 'LaTeX' );

% % behavior of H flat 
% figure(6)
% subplot(1,2,1)
% mesh(x, z, real(H0))
% title('H_{flat}', 'fontname', 'times')
% colorbar;
% 
% subplot(1,2,2)
% mesh(x, z, real(test))
% title('test function','fontname', 'times') 
% colorbar;

% figure(8)
% mesh(x, z, log10( abs( real( H0 - test ) )))
% title('error','fontname', 'times')
% shading flat; 
% colorbar;



