function [J, D_total, tumor_mask, oar_mask, tumor_penalty, oar_penalty] = RadiotherapyCost(theta1, theta2, W_oar)
    % Evaluates cost of beam angles 
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
    tumor_mask = sqrt((X-xc).^2 + (Y-yc).^2) <= 10;
    oar_mask = sqrt((X-xoar).^2 + (Y-yoar).^2) <= 12;
    
    % --- Penalties ---
    % D_total(mask) returns a 1D vector, so only one sum() is needed
    % target is 20 for tumor, max 5 for OAR
    
    tumor_penalty = sum((D_total(tumor_mask) - 20).^2); 
    
    % heavy penalty if OAR dose goes over 5
    oar_penalty = W_oar * sum(max(0, D_total(oar_mask) - 5).^2); 
    
    % tiny penalty to normal tissue just to keep beams from wandering wildly
    normal_penalty = 0.1 * sum(D_total(~tumor_mask & ~oar_mask).^2); 
    
    J = tumor_penalty + oar_penalty + normal_penalty;
    
    
    % disp(['Cost J: ', num2str(J)])
end