% Solving the Hemholtz equation with sound-soft (Dirichlet) boundary
% conditions

% Assuming a flat surface (z = 0)
% Hemholtz eqn: laplacian(u) + k^2 * u = 0, z > 0
% Dirichlet BC: u = 0 on z = 0
% Fundamental solution (in 2D) is Green's Function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% sound-soft BC's require that u_i(x_n,0) (incident field) + u_r(x_n,0) (reflective field) = 0 for n = 1,...,N

k = 1; % wave number 

x0 = 5;  
z0 = 8;

% compute the incident field on z = 0

u_i = 1j / 4 * besselh( 0, 1, k * sqrt( ( x - x0 ).^2 + z0^2 ) );

% compute the MFS system matrix

[ X, Y ] = ndgrid( x );

H = 1j / 4 * besselh( 0, 1, k * sqrt( ( X - Y ).^2 + delta^2 ) );

% weights are obtained from linear system -u_i = Hw

w = H \ -u_i;
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

% compute the reflected field over the mesh

R = sqrt( ( Xmesh(indx) - x(jndx) ).^2 + ( Zmesh(indx) + delta ).^2 );
u_r = 1j / 4 * besselh( 0, 1, k * R ) * w;

% reshape u_r

u_r = reshape( u_r, Nzgrid, Nxgrid );

%% COMPUTE THE EXACT SOLUTION

u_r_exact = - 1j / 4 ...
    * besselh( 0, 1, k * sqrt( ( Xmesh - x0 ).^2 + ( Zmesh + z0 ).^2 ) );

%% plot 

figure(1)
% subplot( 1, 2, 1 );
pcolor( x_grid, z_grid, real( u_r_exact ) );
shading flat;
colorbar;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( '$u_{flat}$: sound-soft', 'Interpreter', 'LaTeX' , 'fontsize', 24);

figure(2)
% subplot( 1, 2, 2 );
pcolor( x_grid, z_grid, real( u_r ) );
shading flat;
colorbar;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
title( '$u_r$ (MFS)', 'Interpreter', 'LaTeX' , 'fontsize', 24);
% sgtitle('Sound-soft', 'Interpreter', 'LaTeX', 'fontsize', 30 );

figure(3)
% subplot( 1, 3, 3 );
pcolor( x_grid, z_grid, log10( abs( u_r - u_r_exact ) ) );
shading flat;
colorbar;
xlabel( '$x$', 'Interpreter', 'LaTeX' );
ylabel( '$z$', 'Interpreter', 'LaTeX' );
hold on
plot( [-L/2 L/2 ], [ 1 1 ], 'o' )
title( 'Absolute Error: sound-soft', 'Interpreter', 'LaTeX' );

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






