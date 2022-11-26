clear all
close all
clc

% 1) Read files solutions

file1 = 'Heat2D-CG-serial-0000032.dat';   % serial 
file2 = 'Heat2D-CG-parallel-0000032.dat'; % parallel

% Open files
fd1 = fopen(file1, 'r');
fd2 = fopen(file2, 'r');

% Read serial run output file
IMAX   = fscanf(fd1, '%d \n', 1);    % Read IMAX
x      = fscanf(fd1, '%f \n', IMAX); % Read x-coords
y      = fscanf(fd1, '%f \n', IMAX); % Read y-coords

% Init variables
T_s = zeros(IMAX,IMAX); 
T_p = zeros(IMAX,IMAX);

% Read serial solution
for j = 1:IMAX
    T_s(:,j) = fscanf(fd1, '%f \n', IMAX);    
end

% Read parallel run output file
IMAX = fscanf(fd2, '%d', 1);       % Read IMAX
x    = fscanf(fd2, '%f \n', IMAX); % Read x-coords
y    = fscanf(fd2, '%f \n', IMAX); % Read y-coords

% Read parallel solution
for j = 1:IMAX
    T_p(:,j) = fscanf(fd2, '%f \n', IMAX);
end

% Boundaries and initialization
TL    = 100.0;          % Temperature on left side
TR    = 50.0;           % Temperature on right side
xD    = 1.0;            % x-coord for discontinuity
kappa = 1.0;            % Heat conduction coefficient
t     = 0.05;           % final time 
T_e = zeros(IMAX,IMAX); % init matrix

% Find the exact solution
for i = 1:IMAX
   for j = 1:IMAX
      T_e(i,j) = 0.5*(TR+TL)...
          + 0.5*erf( (x(j)-xD)/(2.0*sqrt(kappa*t)) )*(TR-TL);
   end
end

% 2) Computational error 

% 2.1) Comparison between exact and serial solution
errTsTe = abs( T_s(1,:)-T_e(1,:) );
% Normalize comp error
errTsTe_norm = errTsTe/norm(errTsTe);

% 2.2) Comparison between exact and parallel solution
errTpTe = abs( T_p(1,:)-T_e(1,:) );
% Normalize comp error
errTpTe_norm = errTpTe/norm(errTpTe);

% 2.3) Comparison between serial and parallel solution
errTsTp = abs( T_p(1,:) - T_s(1,:) );
% Normalized comp error
errTsTp_norm = errTsTp/norm(errTsTp);


% 3) Plot the 2D results ------------------------------------------------%

% 3.1) Serial solution
figure
surf(x,y,T_s);
titleString = '[Serial] Implicit FTCS Scheme for 2-D Heat Equation';
subtitleString = sprintf('TL = 100 °C, TR = 50 °C, IMAX = %d',IMAX);
title(titleString, subtitleString);
xlabel('x');
ylabel('y');
zlabel('Temperature');
colorbar

% 3.2) Parallel solution
figure
surf(x,y,T_p);
titleString = '[Parallel] Implicit FTCS Scheme for 2-D Heat Equation';
subtitleString = sprintf('TL = 100 °C, TR = 50 °C, IMAX = %d',IMAX);
title(titleString, subtitleString);
xlabel('x');
ylabel('y');
zlabel('Temperature');
colorbar

% 3.3) Exact solution
figure
surf(x,y,T_e);
titleString = '[Exact] Implicit FTCS Scheme for 2-D Heat Equation';
subtitleString = sprintf('TL = 100 °C, TR = 50 °C, IMAX = %d',IMAX);
title(titleString, subtitleString);
xlabel('x');
ylabel('y');
zlabel('Temperature');
colorbar

% 4) Plot the 1D solution ----------------------------------------- %

figure
% 4.1) Exact solution
plot(x, T_e(1,:), '--r', 'LineWidth', 1.1);

% 4.2) Serial run solution
hold on 
plot(x, T_s(1,:), '-b',  'LineWidth', 1.1);

% 4.3) Parallel run solution
hold on
plot(x, T_p(1,:), '-.g', 'LineWidth', 1.1);

t = title('Comparison between 1D solutions');
t.FontSize = 16;
legend('Exact', 'Serial Run', 'Parallel Run')
% 5) Plot the computational error --------------------------------- %

% Comparison of computation error of computed solutions

figure
% Norm error between Exact and Serial solution
plot(x, errTsTe_norm, '--or', 'LineWidth', 1.5); 

% Norm error between Parallel and Serial solution
hold on
plot(x, errTpTe_norm, ':*b',  'LineWidth', 1.5);

% Norm error between Parallel and Serial solution
hold on
plot(x, errTsTp_norm, '--xg', 'LineWidth', 1.5);

t = title('Comparison of computational error of computed solutions');
t.FontSize = 14;
legend('Serial vs Exact', 'Parallel vs Serial', 'Serial vs Parallel');

hold off



