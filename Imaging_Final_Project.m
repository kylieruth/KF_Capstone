clear; close all

% Imaging Final Project
% this was before we did PCA

%% FIGURE PARAMETERS

set( 0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0 );

%% MFS

% define the grid 
N = 512;           % must be even 
L = 400;
dx = L / N;
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

% wave number 
k = 2;

% reflectivity 
rho = 10*1j; 

% initial source point
x0 = 30; 
z0 = 30;

% location of target - "point-like"
y0_1 = 10;
y0_2 = 15;

%% COMPUTE INCIDENT FIELD 

A = ( x - x0 ).^2 + ( z - z0 ).^2;

u_i = 1j / 4 * besselh( 0, 1, k * sqrt( A ) );

du_i = -1j * k / 4 * nu_x .* ( x - x0 ) ./ sqrt(A) .* besselh( 1, 1, k * sqrt(A) ) ...
    -1j * k / 4 * nu_z .* ( z - z0 ) ./ sqrt(A) .* besselh( 1, 1, k * sqrt(A) );

%% COMPUTE MFS SYSTEM MATRICES/VECTORS

[ ii, jj ] = ndgrid( (1:N) );    % ii: surface points, jj: source points


% first "block" - reflected field

R = sqrt( ( x(ii) - x(jj) ).^2 + ( z(ii) - z(jj) + dx ).^2 );

% sound-soft
H = 1j / 4 * besselh( 0, 1, k * R );

% sound-hard
H_prime = ( -1j * k / 4  * nu_x(ii) .* (  x(ii) - x(jj)  ) ./ R ...
    .* besselh( 1, 1, k * R ) ) + ( -1j * k / 4 * nu_z(ii) ...
    .* (  z(ii) - z(jj) + dx ) ./ R .* besselh( 1, 1, k * R ) );


% second "block" - scattered field

J = sqrt( ( x - y0_1 ).^2 + ( z - y0_2 ).^2 );

% sound-soft
b = 1j / 4 * besselh( 0, 1, k * J );  

% sound-hard 
db = -1j * k / 4 * nu_x .* ( x - y0_1 ) ./ J .* besselh( 1, 1, k * J ) ...
    -1j * k / 4 * nu_z .* ( z - y0_2 ) ./ J .* besselh( 1, 1, k * J );


% third "block" - self-consistent argument

P = sqrt( ( y0_1 - x ).^2 + ( y0_2 - z + dx ).^2 );

c = 1j / 4 * besselh( 0, 1, k * P );    % sound-soft and sound-hard same; 
                                        % (not evaluated on boundary) 
c = c.';    % want row vector 


% right-hand side %

M = sqrt( ( x0 - y0_1 ).^2 + ( z0 - y0_2 ).^2 );
u_i_y0 = 1j / 4 * besselh( 0, 1, k * M );   % sound-soft = sound-hard

%% SOLVE SYSTEM 

% sound-soft
ss_matrix = [ H  rho * b  ;  c  -1 ];
ss_rhs = [ -u_i ; -u_i_y0 ];

w_ss = ss_matrix \ ss_rhs;
u_e_ss = w_ss(end);                 % exciting field 

% sound-hard
sh_matrix = [ H_prime  rho * db  ;  c  -1 ];
sh_rhs = [ -du_i ; -u_i_y0 ];

w_sh = sh_matrix \ sh_rhs;
u_e_sh = w_sh(end);                 % exciting field 

%% MEASUREMENT MATRIX 

% create grid points for plot 
Nxgrid = 101;
Nzgrid = 201;

x_grid = linspace(-100,100,Nxgrid);
z_grid = linspace(1,100,Nzgrid);

% initialize
D_ss = zeros(Na, Nf);
D_sh = zeros(Na, Nf);

I_ss = zeros(Nzgrid, Nxgrid);
I_sh = zeros(Nzgrid, Nxgrid); 

% mesh grid
[ Xmesh, Zmesh ] = meshgrid( x_grid, z_grid );

% useful index arrays
[ indx, jndx ] = ndgrid( (1:Nxgrid*Nzgrid), (1:N) );


% % calculate measurements 
% R = sqrt( ( Xmesh(indx) - x(jndx) ).^2 + ( Zmesh(indx) - z(jndx) + dx ).^2 );
% J = sqrt( ( Xmesh(:) - y0_1 ).^2 + ( Zmesh(:) - y0_2 ).^2 );
% 
% measurement_ss = 1j / 4 * besselh( 0, 1, k * R ) * w_ss(1:end-1) + ...
% rho * 1j / 4 * besselh( 0, 1, k * J ) * u_e_ss;    % sound-soft
% 
% measurement_sh = 1j / 4 * besselh( 0, 1, k * R ) * w_sh(1:end-1) + ...
% rho * 1j / 4 * besselh( 0, 1, k * J ) * u_e_sh;    % sound-hard
% 
% measurement_ss = reshape( measurement_ss, Nzgrid, Nxgrid );
% measurement_sh = reshape( measurement_sh, Nzgrid, Nxgrid );


% filling measurement matrix 
for m = 1:Nf
    k = k0(m); 

    for n = 1:Na
        x_n = xa(n);
        z_n = za; 

        R = sqrt((x - x_n).^2 + (z - z_n).^2);
        J = sqrt((x_n - y0_1).^2 + (z_n - y0_2).^2);

        % sound-soft
        m_ss = 1j / 4 * besselh(0, 1, k * R).' * w_ss(1:end-1) + ...
            rho * 1j / 4 * besselh(0, 1, k * J) * u_e_ss;

        % sound-hard
        m_sh = 1j / 4 * besselh(0, 1, k * R).' * w_sh(1:end-1) + ...
            rho * 1j / 4 * besselh(0, 1, k * J) * u_e_sh;

        % total scattered field
        D_ss(n, m) = m_ss;
        D_sh(n, m) = m_sh; 

        % compute KM term over grid
        R = sqrt((Xmesh - x_n).^2 + (L - Zmesh).^2);

        I_ss = I_ss + D_ss(n, m) * exp(-1j * 2 * k * R);
        I_sh = I_sh + D_sh(n, m) * exp(-1j * 2 * k * R); 

    end
end

%% PLOTS

% figure(1)
% 
% subplot(1, 2, 1)
% pcolor( x_grid, z_grid, real( measurement_ss ) );
% shading flat;
% colorbar;      % gives height of surface
% hold on
% plot( y0_1, y0_2, 'r+' );
% hold off
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% 
% subplot(1, 2, 2)
% surf( x_grid, z_grid, real( measurement_ss ) );
% shading flat;
% colorbar;      % gives height of surface
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% zlabel( '$u(x,z)$', 'Interpreter', 'Latex');
% sgtitle('Sound-Soft', 'Interpreter', 'Latex');
% 
% figure(2)
% 
% subplot(1, 2, 1)
% pcolor( x_grid, z_grid, real( measurement_sh ) );
% shading flat;
% colorbar;      % gives height of surface
% hold on
% plot( y0_1, y0_2, 'r+' );
% hold off
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% 
% subplot(1, 2, 2)
% surf( x_grid, z_grid, real( measurement_sh ) );
% shading flat;
% colorbar;      % gives height of surface
% xlabel( '$x$', 'Interpreter', 'LaTeX' );
% ylabel( '$z$', 'Interpreter', 'LaTeX' );
% zlabel( '$u(x,z)$', 'Interpreter', 'Latex');
% sgtitle('Sound-Hard', 'Interpreter', 'Latex');

figure(3)
subplot(1, 2, 1);
imagesc(x_grid, z_grid, abs(I_ss));
hold on
plot( y0_1, y0_2, 'r+' );
hold off
title('Sound-Soft', 'Interpreter', 'LaTeX');
xlabel('$x$', 'Interpreter', 'LaTeX'); 
ylabel('$z$', 'Interpreter', 'LaTeX');
colorbar;

subplot(1, 2, 2);
imagesc(x_grid, z_grid, abs(I_sh));
hold on
plot( y0_1, y0_2, 'r+' );
hold off
title('Sound-Hard', 'Interpreter', 'LaTeX');
xlabel('$x$', 'Interpreter', 'LaTeX'); 
ylabel('$z$', 'Interpreter', 'LaTeX');
colorbar;

figure(4)
subplot(1,2,1)
surf( x_grid, z_grid, abs( I_ss ) );
title('Sound-Soft', 'Interpreter', 'LaTeX');
xlabel('$x$', 'Interpreter', 'LaTeX'); 
ylabel('$z$', 'Interpreter', 'LaTeX');
zlabel('$u(x,z)$', 'Interpreter', 'LaTeX'); 
colorbar;
shading flat;

subplot(1,2,2)
surf( x_grid, z_grid, abs( I_sh ) );
title('Sound-Hard', 'Interpreter', 'LaTeX');
xlabel('$x$', 'Interpreter', 'LaTeX'); 
ylabel('$z$', 'Interpreter', 'LaTeX');
zlabel('$u(x,z)$', 'Interpreter', 'LaTeX'); 
colorbar;
shading flat;

