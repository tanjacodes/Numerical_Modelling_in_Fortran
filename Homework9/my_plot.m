% Load data from files
load T_output.dat
load S_output.dat
load W_output.dat
% Load the dynamically generated data
data = load('time_vmax.dat');
time = data(:, 1);  % First column: cumulative time
vmax = data(:, 2);  % Second column: Vmax

% Create a subplot (1 row, 2 columns)
figure;

% First contour plot
subplot(4,1,1);
contourf(T_output);
title('T');
colorbar;
subplot(4,1,2);
contourf(W_output);
title('W');
colorbar;
subplot(4,1,3);
contourf(S_output);
title('S');
colorbar;
subplot(4,1,4);
% Plot the results
plot(time, vmax, 'LineWidth', 1.5);
xlabel('Time');            % Label for x-axis
ylabel('V_{max}');         % Label for y-axis
title('V_{max} vs Time');  % Title