clear all
close all
clc

% Open the file in read mode
fd = fopen('Heat2D-0002202.dat', 'r');

% Read IMAX
IMAX = fscanf(fd,'%f \n',1);

% Variables initialization
T(IMAX,IMAX) = 0.0;
x(IMAX)      = 0.0;
y(IMAX)      = 0.0;

% Read x-coords
x = fscanf(fd, '%f \n', IMAX);

% Read y-coords
y = fscanf(fd, '%f \n', IMAX);

% Read the temperature
for j = 1:IMAX 
    T(:,j) = fscanf(fd, '%f \n', IMAX);
end

% Close the file
fclose(fd);

% Ploit the result and make the plot prettier 
surf(x,y,T);
titleString = 'FTCS Scheme for 2-D Heat Equation';
stitleString = sprintf('TL = 100 °C, TR = 50 °C; N = %d', IMAX);
[t,s] = title(titleString, stitleString);
t.FontSize = 16;
s.FontSize = 11;
s.FontAngle = 'italic';
xlabel('x')
ylabel('y')
zlabel('Temperature')
colorbar



