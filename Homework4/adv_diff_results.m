% Load data from files
load f_vcycle_random.dat
load u_vcycle_random.dat

% Create a subplot (1 row, 2 columns)
figure;

% First contour plot
subplot(1, 2, 1);
contourf(u_vcycle_random);
title('u vcycle random');
colorbar;

% Second contour plot
subplot(1, 2, 2);
contourf(adv_diff_result);
title('Adv Diff Result');
colorbar;