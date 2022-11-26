%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Exact solution of the 2-D heat conduction equation.             %
%                Written by Enrico Bregoli, 08-11-2022                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%
% Heat conduction coefficient
kappa = 1.0;
%
% Set boundary conditions
TL = 100.0;     % left temperature
TR = 50.0;      % right temperature
%
% Set final time
t = 0.05;
%
% Set computational domain 
xL   = 0.0;     % left boundary
xR   = 2.0;     % right boundary
yB   = 0.0;     % bottom boudary
yT   = 2.0;     % top boundary
xD   = 1.0;     % x-coord of discontinuity
IMAX = 50;      % total number of cells
x = linspace(xL,xR,IMAX); % x-coordinates of cells barycenters
y = linspace(yB,yT,IMAX); % y-coordinates of cells barycenters
%
Te = zeros(IMAX,IMAX); % init matrix
% Exact solution
for i = 1:IMAX
   for j = 1:IMAX
      Te(i,j) = 0.5*(TR+TL)...
          + 0.5*erf( (x(j)-xD)/(2.0*sqrt(kappa*t)) )*(TR-TL);
   end
end
%
% Plot the results and add some stuff to make plot looks prettier :-)
surf(x,y,Te)
% title and subtitle
[t,s] = title( 'Exact solution for 2-D Heat Equation', ...
                sprintf('TL = 100 °C, TR = 50 °C; N = %d',IMAX) );
t.FontSize = 16;
s.FontSize = 10;
s.FontAngle = 'italic';
% add colorbar
colorbar
% add labels to axis
xlabel('x')
ylabel('y')
zlabel('Temperature')


