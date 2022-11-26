clear
close all
clc

% Plot computational time -------------------------------------------------

% Open files 
% Explicit solver results
fdes = fopen("./Heat2Dex_ExecTime/AvgExecTime_Heat2D_EX_0001.dat", "r"); 
% Implicit solver results
fdis = fopen("./Heat2Dcg_ExecTime/AvgExecTIme_Heat2D_CG_0001.dat", "r");
% Explicit solver in parallel results 
fdep = fopen("./Heat2Dex_ExecTime/AvgExecTime_Heat2D_EX_0004.dat", "r");
% Implicit solver in parallel results
fdip = fopen("./Heat2Dcg_ExecTime/AvgExecTIme_Heat2D_CG_0004.dat", "r");

res_es = fscanf(fdes, '%f', [2,inf])'; % read explicit serial results
res_is = fscanf(fdis, '%f', [2,inf])'; % read implicit serial results
res_ep = fscanf(fdep, "%f", [2,inf])'; % read explicit parallel results
res_ip = fscanf(fdip, "%f", [2,inf])'; % read implicit parallel results

% close files
fclose(fdes);
fclose(fdis);
fclose(fdep);
fclose(fdip);

x  = res_es(:, 1); % Get x-coords (same for serial and parallel)
x = x .* x;
yes = res_es(:, 2); % Get y-coords for explicit serial
yis = res_is(:, 2); % Get y-coords for implicit serial
yep = res_ep(:, 2); % Get y-coords for explicit parallel
yip = res_ip(:, 2); % Get y-coords for implicit parallel

% Plot the explicit serial results, ...
figure
plot(x, yes, 'r', 'LineWidth', 1.5);
% axis label
xlabel('Dimensione della griglia (numero di elementi)');
ylabel('Tempo di computazione medio (secondi)');
% title
[t, s] = title('Tempo medio di computazione [Heat 2-D]', 'Intel(R) Core(TM) i7-8750H');
t.FontSize = 14;

% ... implicit serial results, ...
hold on
plot(x, yis, 'b', 'LineWidth', 1.5);

% ... explicit parallel results, ...
hold on
plot(x, yep, 'g', 'LineWidth', 1.5);

% implicit serial results, ...
hold on
plot(x, yip, 'm', "LineWidth", 1.5);

% axix size
axis([min(x) max(x) min(yes) max(yes)]);
% legend
legend('Risolutore esplicito (1 CPU)', 'Risolutore implicito C.G. (1 CPU)', ...
       'Risolutore esplicito (4 CPU)', 'Risolutore implicito C.G. (4 CPU)');
hold off
