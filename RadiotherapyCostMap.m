function Z = RadiotherapyCostMap(THETA1, THETA2, W_oar)
    % helper function to generate the 3D surface plot
    % runs the cost function over the entire grid
    
    Z = zeros(size(THETA1));
    
    % nested loop is slow but it's just for visualization so it's fine
    for i = 1:size(THETA1, 1)
        for j = 1:size(THETA1, 2)
            Z(i,j) = RadiotherapyCost(THETA1(i,j), THETA2(i,j), W_oar);
        end
    end
end