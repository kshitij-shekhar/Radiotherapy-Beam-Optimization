%% BME Optimization Methods Project
% Author: Kshitij Shekhar, Dhanush Kumar (Adapted from Prof. E. Schkommodau)
% Description: Hybrid Algorithm (PSO + Nelder-Mead) for Radiotherapy Angle Opt.

clear; close all; clc;


W_oar = 20; % penalty weight for OAR 
bounds = [0, pi]; % keep beams between 0 and 180 degrees

f_cost = @(t1, t2) RadiotherapyCost(t1, t2, W_oar);

%% Phase 1: PSO Global Search
disp('Running Phase 1: Swarm Search...');

% PSO params (don't need a massive swarm for just 2 variables)
pop = 5; 
c1 = 0.6; c2 = 0.6; c3 = 0.3; %weight for each particle's best known pos, pull towards swarm's best known pos, 
max_iter_pso = 6;

[g_best, pso_cost] = run_PSO(f_cost, bounds, pop, c1, c2, c3, max_iter_pso);
% g_best is PSO optimized [theta1,theta2], pso_cost is scalar best cost of
% PSO optimized [theta1,theta2]

fprintf('PSO rough guess: theta1 = %.2f, theta2 = %.2f (Cost: %.2f)\n', g_best(1), g_best(2), pso_cost);
disp('Press any key to start Nelder-Mead tumble...');
pause; 

%% Phase 2: Nelder-Mead Refinement
disp('Starting Phase 2: Downhill Simplex...');

% build the initial triangle  around the PSO's best guess
step = 0.06;%0.2; 
triangle = [g_best; 
            g_best + [step, 0]; 
            g_best + [0, step]];

% initial triangle coordinates : Point 1 - [theta1,theta2], Point 2 :
% [theta1+0.2,theta2], Point 3 : [theta1,theta2+0.2] 
% triangle coordinates are in rad

%Extension : clamp Nelder Mead simplex : 
triangle = max(min(triangle, bounds(2)), bounds(1));
% ----


alpha = 1; beta = 2; gamma = 0.5; delta = 0.5;
eps_break = 1e-4;
max_iter_nm = 20; 

% setup the 3D surface plot for visualization
[x1, x2] = meshgrid(0:0.01:pi, 0:0.01:pi);
% [x1, x2] = meshgrid(0:0.1:pi, 0:0.1:pi);
f_map = RadiotherapyCostMap(x1, x2, W_oar); 

%Extension - adding contour plot to view NM in action
fig2D = figure;
contourf(x1, x2, f_map, 30, 'LineStyle', 'none');
hold on;
contour(x1, x2, f_map, 20, 'k');
colormap parula;
colorbar;
xlabel('Beam Angle 1');
ylabel('Beam Angle 2');
title('Nelder-Mead on Cost Contours');
axis tight;

fig3D = figure;
hold on;
surf(x1, x2, f_map, 'EdgeColor', 'none');
colormap parula; view(3);
xlabel('Beam Angle 1'); ylabel('Beam Angle 2'); zlabel('Penalty');

% ----



% figure;
% contourf(x1, x2, f_map, 30, 'LineStyle', 'none');
% hold on;
% contour(x1, x2, f_map, 20, 'k');
% colormap parula;
% colorbar;
% xlabel('Beam Angle 1');
% ylabel('Beam Angle 2');
% title('Nelder-Mead on Cost Contours');
% axis tight;
% 
% figure; hold on;
% surf(x1, x2, f_map, 'EdgeColor', 'none'); 
% colormap parula; view(3);
% xlabel('Beam Angle 1'); ylabel('Beam Angle 2'); zlabel('Penalty');


iteration = 1;
epsilon = 1000;

while (iteration < max_iter_nm) && (epsilon > eps_break)
    % evaluate the 3 corners
    t_vals = [f_cost(triangle(1,1), triangle(1,2)); ... % evaluate cost for first point of triangle : it's theta1 and theta2 vals
              f_cost(triangle(2,1), triangle(2,2)); ...
              f_cost(triangle(3,1), triangle(3,2))];
    %t_vals is a 3x1 matrix of scalar values

    % sort points: best, intermediate, worst 
    [~, sorted_idx] = sort(t_vals);
    x_b = triangle(sorted_idx(1), :); 
    x_i = triangle(sorted_idx(2), :); 
    x_w = triangle(sorted_idx(3), :); 

    % plot current triangle
    c_line = [rem(iteration,2), iteration/max_iter_nm, 0];

    %Extension : contour overlay
    % 2D contour overlay of current simplex
    figure(fig2D);
    plot([x_b(1) x_i(1) x_w(1) x_b(1)], ...
         [x_b(2) x_i(2) x_w(2) x_b(2)], ...
         '-', 'Color', c_line, 'LineWidth', 2);

    plot(x_b(1), x_b(2), 'og', 'MarkerFaceColor', 'g', 'MarkerSize', 6);

    figure(fig3D);
    % -----
    line([x_b(1),x_i(1)], [x_b(2),x_i(2)], [f_cost(x_b(1),x_b(2)), f_cost(x_i(1),x_i(2))], 'LineWidth',2,'Color', c_line);
    line([x_i(1),x_w(1)], [x_i(2),x_w(2)], [f_cost(x_i(1),x_i(2)), f_cost(x_w(1),x_w(2))], 'LineWidth',2,'Color', c_line);
    line([x_w(1),x_b(1)], [x_w(2),x_b(2)], [f_cost(x_w(1),x_w(2)), f_cost(x_b(1),x_b(2))], 'LineWidth',2,'Color', c_line);
    plot3(x_b(1), x_b(2), f_cost(x_b(1),x_b(2)), 'og', 'LineWidth', 2);

    title(sprintf('Nelder-Mead Iteration %d - Press any key', iteration));
    pause; 

    % geometric tumble logic
    x_mean = 0.5 * (x_i + x_b); % midpoint of the edge opposite the worst corner
    x_r = x_mean + alpha*(x_mean - x_w);  
    % Extension:
    x_r = max(min(x_r, bounds(2)), bounds(1));

    if f_cost(x_b(1),x_b(2)) <= f_cost(x_r(1),x_r(2)) && f_cost(x_r(1),x_r(2)) < f_cost(x_i(1),x_i(2)) 
         x_w = x_r; % reflect
    elseif f_cost(x_r(1),x_r(2)) < f_cost(x_b(1),x_b(2))
        
        x_e = x_mean + beta * (x_r - x_mean); %expand and take an even bigger step in case the reflected point is better than best point
        x_e = max(min(x_e, bounds(2)), bounds(1));

        if f_cost(x_e(1),x_e(2)) < f_cost(x_r(1),x_r(2)) 
            x_w = x_e; % expand
        else 
            x_w = x_r;
        end
    elseif f_cost(x_i(1),x_i(2)) <= f_cost(x_r(1),x_r(2)) && f_cost(x_r(1),x_r(2)) < f_cost(x_w(1),x_w(2))   
        x_oc = x_mean + gamma * (x_r - x_mean);
        x_oc = max(min(x_oc, bounds(2)), bounds(1));

        if f_cost(x_oc(1),x_oc(2)) <= f_cost(x_r(1),x_r(2))
            x_w = x_oc; % outer contract
        else
            x_i = x_b + delta * (x_i - x_b); x_w = x_b + delta * (x_w - x_b); % shrink
        end
    elseif f_cost(x_r(1),x_r(2)) >= f_cost(x_w(1),x_w(2))
        x_ic = x_mean - gamma * (x_r - x_mean);
        x_ic = max(min(x_ic, bounds(2)), bounds(1));

        if f_cost(x_ic(1),x_ic(2)) < f_cost(x_r(1),x_r(2))
            x_w = x_ic; % inner contract
        else
            x_i = x_b + delta * (x_i - x_b); x_w = x_b + delta * (x_w - x_b);
        end
    end

    epsilon = norm(abs(x_w - x_i) + abs(x_w - x_b) + abs(x_b - x_i));
    triangle = [x_b; x_i; x_w];
    iteration = iteration + 1;
end 
nm_cost = f_cost(x_b(1), x_b(2));
fprintf('PSO cost: %.6f\n', pso_cost);
fprintf('NM final cost: %.6f\n', nm_cost);
fprintf('NM improvement: %.6f\n', pso_cost - nm_cost);
fprintf('Percent improvement: %.4f%%\n', 100*(pso_cost - nm_cost)/pso_cost);

plot3(x_b(1), x_b(2), f_cost(x_b(1),x_b(2)), 'ow', 'MarkerSize', 8, 'LineWidth', 2);
title('Phase 2 Complete');
fprintf('Final optimized angles: theta1 = %.3f, theta2 = %.3f\n', x_b(1), x_b(2));

%% Clinical Visuals (Heatmap & DVH)
% actual dose map for the final angles
[~, dose_map, t_mask, o_mask] = RadiotherapyCost(x_b(1), x_b(2), W_oar);

figure;
contourf(dose_map, 30, 'LineStyle', 'none'); colormap jet; hold on;
contour(t_mask, [1 1], 'r', 'LineWidth', 2); % Tumor
contour(o_mask, [1 1], 'w', 'LineWidth', 2); % OAR
title('Final Radiation Dose Distribution'); 

figure; hold on;
t_doses = sort(dose_map(t_mask), 'descend');
o_doses = sort(dose_map(o_mask), 'descend');
plot(t_doses, linspace(100, 0, length(t_doses)), 'r', 'LineWidth', 2); 
plot(o_doses, linspace(100, 0, length(o_doses)), 'b', 'LineWidth', 2);     
xlabel('Dose'); ylabel('Volume (%)');
title('Dose-Volume Histogram'); legend('Tumor', 'OAR'); grid on;

%% Pareto Front (Sensitivity Analysis)
disp('Running Pareto Sweep ');
weights = [1, 5, 20, 50, 150];
t_errs = zeros(1, length(weights));
o_errs = zeros(1, length(weights));

for i = 1:length(weights)
    fw = @(t1, t2) RadiotherapyCost(t1, t2, weights(i));

    % silent PSO
    [rough_best, ~] = run_PSO(fw, bounds, pop, c1, c2, c3, 10);

    % using built-in fminsearch here so the loop doesn't take 10 years

    exact_best = fminsearch(@(x) fw(x(1), x(2)), rough_best);

    [~, ~, ~, ~, t_pen, o_pen] = RadiotherapyCost(exact_best(1), exact_best(2), weights(i));
    t_errs(i) = t_pen;
    o_errs(i) = o_pen / weights(i); % strip the weight multiplier to get raw dose error
end

figure;
plot(t_errs, o_errs, '-ro', 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
xlabel('Tumor Underdose Error'); ylabel('OAR Overdose Error');
title('Pareto Front: W_{oar} Sensitivity'); grid on;

% %% BME Optimization Methods Project
% % Description: Hybrid Algorithm for Radiotherapy Angle Opt. (2 OARs)
% 
% clear; close all; clc;
% 
% % We now define two separate weights
% W_oar1 = 50; 
% W_oar2 = 80; % We made OAR 2 even more critical
% bounds = [0, pi]; 
% 
% % Pass both weights into the anonymous function
% f_cost = @(t1, t2) RadiotherapyCost(t1, t2, W_oar1, W_oar2);
% 
% %% Phase 1: PSO Global Search
% disp('Running Phase 1: Swarm Search...');
% pop = 5; 
% c1 = 0.6; c2 = 0.6; c3 = 0.3; 
% max_iter_pso = 15;
% 
% [g_best, pso_cost] = run_PSO(f_cost, bounds, pop, c1, c2, c3, max_iter_pso);
% 
% fprintf('PSO rough guess: theta1 = %.2f, theta2 = %.2f (Cost: %.2f)\n', g_best(1), g_best(2), pso_cost);
% disp('Press any key to start Nelder-Mead tumble...');
% pause; 
% 
% %% Phase 2: Nelder-Mead Refinement
% disp('Starting Phase 2: Downhill Simplex...');
% 
% % fminsearch uses the exact Nelder-Mead mathematics
% options = optimset('Display', 'iter', 'TolX', 1e-4);
% [best_angles, final_cost] = fminsearch(@(x) f_cost(x(1), x(2)), g_best, options);
% 
% fprintf('\nOPTIMIZATION COMPLETE\n');
% fprintf('Final Theta 1: %.2f degrees\n', rad2deg(best_angles(1)));
% fprintf('Final Theta 2: %.2f degrees\n', rad2deg(best_angles(2)));
% 
% %% Clinical Visualizations
% % Extract final map and masks (note the extra output for OAR2)
% [~, dose_map, t_mask, o1_mask, o2_mask] = RadiotherapyCost(best_angles(1), best_angles(2), W_oar1, W_oar2);
% 
% % 1. Heatmap
% figure; 
% contourf(dose_map, 30, 'LineStyle', 'none'); colormap jet; hold on;
% contour(t_mask, [1 1], 'r', 'LineWidth', 2); % Tumor (Red)
% contour(o1_mask, [1 1], 'w', 'LineWidth', 2); % OAR 1 (White)
% contour(o2_mask, [1 1], 'm', 'LineWidth', 2); % OAR 2 (Magenta)
% title('Final Radiation Dose Distribution'); 
% 
% % 2. Dose-Volume Histogram (DVH)
% figure; hold on;
% t_doses = sort(dose_map(t_mask), 'descend');
% o1_doses = sort(dose_map(o1_mask), 'descend');
% o2_doses = sort(dose_map(o2_mask), 'descend');
% 
% plot(t_doses, linspace(100, 0, length(t_doses)), 'r', 'LineWidth', 2); 
% plot(o1_doses, linspace(100, 0, length(o1_doses)), 'b', 'LineWidth', 2);     
% plot(o2_doses, linspace(100, 0, length(o2_doses)), 'm', 'LineWidth', 2);  % New line for OAR 2   
% 
% xlabel('Dose'); ylabel('Volume (%)');
% title('Dose-Volume Histogram'); 
% legend('Tumor', 'OAR 1 (Spinal Cord)', 'OAR 2 (Heart)'); grid on;
% 
% %% Pareto Front (Sensitivity Analysis)
% disp('Running Pareto Sweep ');
% weights = [1, 5, 20, 50, 150];
% t_errs = zeros(1, length(weights));
% total_oar_errs = zeros(1, length(weights)); % We now sum the OAR errors to keep the graph 2D
% 
% for i = 1:length(weights)
%     % Sweep both weights simultaneously
%     fw = @(t1, t2) RadiotherapyCost(t1, t2, weights(i), weights(i));
%     [rough_best, ~] = run_PSO(fw, bounds, pop, c1, c2, c3, 10);
%     best_ang = fminsearch(@(x) fw(x(1), x(2)), rough_best);
% 
%     [~, ~, ~, ~, ~, t_pen, o1_pen, o2_pen] = RadiotherapyCost(best_ang(1), best_ang(2), weights(i), weights(i));
% 
%     t_errs(i) = t_pen;
%     % Combine OAR 1 and OAR 2 penalties for the Y-axis
%     total_oar_errs(i) = o1_pen + o2_pen; 
% end
% 
% figure;
% plot(t_errs, total_oar_errs, '-bo', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% xlabel('Tumor Error (Underdose)'); ylabel('Combined OAR Error (Toxicity)');
% title('Pareto Trade-Off Front'); grid on;