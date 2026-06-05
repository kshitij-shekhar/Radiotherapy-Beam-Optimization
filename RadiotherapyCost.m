function [J, D_total, tumor_mask, oar_mask, tumor_penalty, oar_penalty] = RadiotherapyCost(theta1, theta2, W_oar)
    % Evaluates cost of beam angles for both PSO and NM
    % theta1, theta2: input angles (rad)
    % W_oar: penalty weight for organ at risk

    % setup grid (100x100 tissue slice)
    [X, Y] = meshgrid(1:100, 1:100);

    % tissue centers
    xc = 50; yc = 50;     % tumor center
    xoar = 35; yoar = 35; % OAR

    % beam params 
    I_max = 10;     % max dose at center ray
    sigma = 5;      % beam spread
    denom = 2 * sigma^2;

    % --- Beam 1 ---
    % perp distance to beam line
    d1 = abs((X - xc)*sin(theta1) - (Y - yc)*cos(theta1));
    D1 = I_max * exp(-(d1.^2) / denom); 

    % --- Beam 2 ---
    d2 = abs((X - xc)*sin(theta2) - (Y - yc)*cos(theta2));
    D2 = I_max * exp(-(d2.^2) / denom);

    % superposition
    D_total = D1 + D2; 

    % generate logical masks for the regions
    tumor_mask = sqrt((X-xc).^2 + (Y-yc).^2) <= 10; %defining the tumor regions in the tissue grid
    oar_mask = sqrt((X-xoar).^2 + (Y-yoar).^2) <= 12; %defining the oar regions in the tissue grid

    % --- Penalties ---
    % D_total(mask) returns a 1D vector, so only one sum() is needed
    % target is 20 for tumor, max 5 for OAR

    tumor_penalty = sum((D_total(tumor_mask) - 20).^2); % D_total(tumor_mask) gives all dose values inside the tumor region

    % heavy penalty if OAR dose goes over 5
    oar_penalty = W_oar * sum(max(0, D_total(oar_mask) - 5).^2); 

    %------------
    %Extension : adding 1 more OAR
    % xoar2 = 65; yoar2 = 40;
    % oar2_mask = sqrt((X-xoar2).^2 + (Y-yoar2).^2) <= 10;
    % 
    % oar_penalty2 = 0.8*W_oar * sum(max(0, D_total(oar2_mask) - 5).^2);
    %------------

    % penalty to prevent beams from irradiating normal tissue unecessarily
    normal_penalty = 0.1 * sum(D_total(~tumor_mask & ~oar_mask).^2); 

    J = tumor_penalty + oar_penalty  + normal_penalty; %+ oar_penalty2;


    % disp(['Cost J: ', num2str(J)])
end

% function [J, D_total, tumor_mask, oar1_mask, oar2_mask, tumor_penalty, oar1_penalty, oar2_penalty] = RadiotherapyCost(theta1, theta2, W_oar1, W_oar2)
%     % Evaluates cost of beam angles with TWO Organs at Risk
% 
%     % setup grid (100x100 tissue slice)
%     [X, Y] = meshgrid(1:100, 1:100);
% 
%     % tissue centers
%     xc = 50; yc = 50;       % Tumor center
%     xoar1 = 35; yoar1 = 35; % OAR 1 (e.g., Spinal Cord)
%     xoar2 = 70; yoar2 = 65; % OAR 2 (e.g., Heart) - Placed on the opposite side
% 
%     % beam params 
%     I_max = 10;     
%     sigma = 3;      
%     denom = 2 * sigma^2;
% 
%     % --- Beam 1 ---
%     d1 = abs((X - xc)*sin(theta1) - (Y - yc)*cos(theta1));
%     D1 = I_max * exp(-(d1.^2) / denom); 
% 
%     % --- Beam 2 ---
%     d2 = abs((X - xc)*sin(theta2) - (Y - yc)*cos(theta2));
%     D2 = I_max * exp(-(d2.^2) / denom);
% 
%     % superposition
%     D_total = D1 + D2; 
% 
%     % generate logical masks
%     tumor_mask = sqrt((X-xc).^2 + (Y-yc).^2) <= 10; 
%     oar1_mask = sqrt((X-xoar1).^2 + (Y-yoar1).^2) <= 12; 
%     oar2_mask = sqrt((X-xoar2).^2 + (Y-yoar2).^2) <= 9; % OAR 2 is slightly smaller
% 
%     % --- Penalties ---
%     tumor_penalty = sum((D_total(tumor_mask) - 20).^2); 
% 
%     % heavy penalty if OAR 1 dose goes over 5
%     oar1_penalty = W_oar1 * sum(max(0, D_total(oar1_mask) - 5).^2); 
% 
%     % heavy penalty if OAR 2 dose goes over 4 (Stricter limit!)
%     oar2_penalty = W_oar2 * sum(max(0, D_total(oar2_mask) - 4).^2);
% 
%     % Normal tissue sweep (ignores tumor and BOTH OARs)
%     n_pen = 0.1 * sum(sum(D_total(~tumor_mask & ~oar1_mask & ~oar2_mask).^2));
% 
%     % Total Cost
%     J = tumor_penalty + oar1_penalty + oar2_penalty + n_pen;
% end