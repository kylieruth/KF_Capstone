% Solving the Hemholtz equation with sound-hard (Neumann) boundary
% conditions

clear;

%% FIGURE PARAMETERS

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0);

%% Method of Fund. Solutions (MFS)

% define the grid 

N = 100;  
L = 50;
delta = 0.1;  

z = -delta;
x = zeros(N,1);

for n = 1:N
    x(n) = (-L / 2) + ((n - 1) * L) / (N - 1);
end

% sound-hard BC's require that dz( u_i(x_n,z) + u_r(x_n,z) ) = 0 for n = 1,...,N

k = 1; % wave constant 

x0 = 5; 
z0 = 8;

% compute the incident field and normal deriviative 

A = ( x - x0 ).^2 + z0^2;
u_i = 1j / 4 * besselh( 0, 1, k * sqrt( A ) );

du_i_dz = 1j * k / 4.0 * z0 ./ sqrt( A ) .* besselh( 1, 1, k * sqrt( A) );

% compute the MFS system matrix

[ X, Y ] = ndgrid( x );

H = 1j / 4 * besselh( 0, 1, k * sqrt( ( X - Y ).^2 + delta^2 ) );

R    = sqrt( ( X - Y ).^2 + delta^2 );
H_dz = - 1j * k / 4.0 * delta ./ R .* besselh( 1, 1, k * R );

% weights are obtained from linear system -u_i = Hw
w = H_dz \ - du_i_dz;
disp(w)

% create grid points for plot 

Nxgrid = 101;
Nzgrid = 201;

x_grid = linspace(-50,50,Nxgrid);
z_grid = linspace(1,50,Nzgrid);

% mesh grid

[ Xmesh, Zmesh ] = meshgrid( x_grid, z_grid );

% useful index arrays

[ indx, jndx ] = ndgrid( (1:Nxgrid*Nzgrid), (1:N) );

% compute the reflected field over the mesh and normal derivative

R = sqrt( ( Xmesh(indx) - x(jndx) ).^2 + ( Zmesh(indx) + delta ).^2 );
u_r = 1j / 4 * besselh( 0, 1, k * R ) * w; 

% reshape

u_r = reshape( u_r, Nzgrid, Nxgrid );

%% COMPUTE THE EXACT SOLUTION

u_r_exact = 1j / 4 ...
    * besselh( 0, 1, k * sqrt( ( Xmesh - x0 ).^2 + ( Zmesh + z0 ).^2 ) );  % made sign changes

%% plot 

figure(1)
% subplot( 1, 2, 1 );
pcolor( x_grid, z_grid, real( u_r_exact ) );
shading flat;
colorbar;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( '$u_{flat}$: sound-hard', 'Interpreter', 'LaTeX' , 'fontsize', 24);

figure(2)
% subplot( 1, 2, 2 );
pcolor( x_grid, z_grid, real( u_r ) );
shading flat;
colorbar;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( '$u_r$ (MFS)', 'Interpreter', 'LaTeX' , 'fontsize', 24);
% sgtitle('Sound-hard', 'Interpreter', 'LaTeX', 'fontsize', 24);

figure(3)
% subplot( 1, 3, 3 );
pcolor( x_grid, z_grid, log10( abs( u_r - u_r_exact ) ) );
shading flat;
colorbar;
hold on;
plot( [-L/2 L/2 ], [ 1 1 ], 'o' )
hold off;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( 'Absolute Error: sound-hard', 'Interpreter', 'LaTeX' );

% figure(2)
% 
% subplot(1,2,1)
% pcolor( x_grid, z_grid, real(u_r) )
% shading flat;
% colorbar;
% 
% subplot(1,2,2)
% surf( x_grid, z_grid, real(u_r) )
% shading flat;
% colorbar;
% 
% figure(3)
% subplot(1,2,1)
% imagesc( abs(H) );
% colorbar;
% 
% subplot(1,2,2)
% imagesc( abs(H_dz) );
% colorbar;