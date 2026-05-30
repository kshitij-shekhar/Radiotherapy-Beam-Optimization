function [g_best_pos, g_best_cost] = run_PSO(cost_func, bounds, pop_size, c1, c2, c3, max_iter)
    % basic PSO to find the general valley before Nelder-Mead takes over
    
    lb = bounds(1); 
    ub = bounds(2); 
    
    % drop particles randomly within the [0, pi] bounds
    curr_pos = lb + (ub - lb) * rand(pop_size, 2);
    vel = zeros(pop_size, 2);
    
    best_local_costs = zeros(pop_size, 1);
    
    % initial eval
    for i = 1:pop_size
        best_local_costs(i) = cost_func(curr_pos(i,1), curr_pos(i,2));
    end
    
    best_local_pos = curr_pos;
    
    % find initial global best
    [g_best_cost, idx] = min(best_local_costs);
    g_best_pos = curr_pos(idx, :);
    
    % main swarm loop
    for iter = 1:max_iter
        r1 = rand(pop_size, 2);
        r2 = rand(pop_size, 2);
        
        % standard PSO velocity update
        vel = c3*vel + c1*r1.*(best_local_pos - curr_pos) + c2*r2.*(repmat(g_best_pos, pop_size, 1) - curr_pos);
                   
        % move particles and clamp them so the beams don't try to spin 360 degrees
        new_pos = curr_pos + vel;
        new_pos = max(min(new_pos, ub), lb);
        
        % evaluate new spots
        for i = 1:pop_size
            new_cost = cost_func(new_pos(i,1), new_pos(i,2));
            
            if new_cost < best_local_costs(i)
                best_local_costs(i) = new_cost;
                best_local_pos(i,:) = new_pos(i,:);
            end
        end
        
        % update swarm's overall best
        [current_min, min_idx] = min(best_local_costs);
        if current_min < g_best_cost
            g_best_cost = current_min;
            g_best_pos = best_local_pos(min_idx, :);
        end
        
        curr_pos = new_pos;
    end
end