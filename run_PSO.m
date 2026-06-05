function [g_best_pos, g_best_cost] = run_PSO(cost_func, bounds, pop_size, c1, c2, c3, max_iter)
    % basic PSO to find the general valley before Nelder-Mead takes over
    
    lb = bounds(1); %lower bound = 0
    ub = bounds(2); %upper bound = pi 
    
    % drop particles randomly within the [0, pi] bounds
    curr_pos = lb + (ub - lb) * rand(pop_size, 2); % rand(pop_size,2) -> a pop_sizex2 matrix of random numbers between 0 and 1
    % lb+(ub-lb) rescales the matrix so that the random numbers lie
    % between ub and lb

    %curr_pos is thus a pop_sizex2 matrix with each row a particle and each
    %column a beam angle for that particular particle

    vel = zeros(pop_size, 2); % initialises each particle's step size and step direction 
    
    best_local_costs = zeros(pop_size, 1); %initialize a column vector with each particle's best cost
    
    % initial eval
    for i = 1:pop_size
        best_local_costs(i) = cost_func(curr_pos(i,1), curr_pos(i,2)); %curr_pos(i,1) = particle i's theta 1
    end
    
    best_local_pos = curr_pos; % stores the current pop_sizex2 matrix of beam angles for each particle as their best
    
    % find initial global best
    [g_best_cost, idx] = min(best_local_costs); %find lowest cost amongst the particles' personal best costs
    g_best_pos = curr_pos(idx, :); %stores lowest cost particle's pos
    
    % ----
    % At this point we have found where every particle currently is, the
    % best place each particle has found, and the best place the whole
    % swarm has found
    % ----
    % main swarm loop
    for iter = 1:max_iter
        r1 = rand(pop_size, 2); %pop_sizex2 random matrix mimicking particle pos matrix gives PSO some stochasticity
        r2 = rand(pop_size, 2);
        
        % standard PSO velocity update
        vel = c3*vel + c1*r1.*(best_local_pos - curr_pos) + c2*r2.*(repmat(g_best_pos, pop_size, 1) - curr_pos);
        % c1*r1 - cognitive term (pulls each particle towards its own best known location)
        % c2*r2 - social term (pulls each particle towards swarm's best
        % known location) 
        % c3*vel - inertia term for stability in particle movement

        % move particles and clamp them so the beams don't try to spin 360 degrees
        new_pos = curr_pos + vel; %new_pos(i) = (theta1+vel1, theta2+vel2)
        new_pos = max(min(new_pos, ub), lb);
        
        % evaluate new spots
        for i = 1:pop_size
            new_cost = cost_func(new_pos(i,1), new_pos(i,2)); % new cost for each particle i
            
            if new_cost < best_local_costs(i)
                best_local_costs(i) = new_cost; %update best local cost of particle i if the new cost is lower than the best
                best_local_pos(i,:) = new_pos(i,:); %update best local pos of the particle as well
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