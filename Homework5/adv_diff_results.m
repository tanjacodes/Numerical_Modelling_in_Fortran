% Load data from files
load T_output.dat
%load u_vcycle_random.dat

% Create a subplot (1 row, 2 columns)
figure;

% First contour plot
%subplot(1, 2, 1);
contourf(T_output);
title('T');
colorbar;

% Second contour plot
%subplot(1, 2, 2);
%contourf(f_vcycle_random);
title('f vcycle random');
colorbar;