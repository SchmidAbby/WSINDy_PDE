function displacement_data = integrate_data(velocity_data, time)
    
    [x_pts, t_pts] = size(velocity_data);

    % initialize matrix for displacement data
    displacement_data = zeros(x_pts, t_pts);

    for i=1:x_pts
        displacement_data(i, :) = cumtrapz(velocity_data(i,:), time);
    end    

end