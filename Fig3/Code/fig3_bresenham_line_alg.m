warning('off')
%%
image_num = 21600;
x_size = 151;
y_size = 101;


%% filter
fil=zeros(26,26);
for l=1:26
    for m=1:26
        if (l-12.5)^2 +(m-11)^2 > 6^2 && (l-12.5)^2 +(m-11)^2 < 9.5^2
            fil(l,m)=1;
        end
    end
end


%%          
pixel_mat_x = [17, 18, 19, 20, 20, 20, 20, 20, 20, 20, 21, 20, 19, 18, 17, 17, 16, 15, 14, 13, 12, 10, 9, 8, 7, 6, 4, 3, 4, 5, 6, 6, 6, 6, 6, 7, 8, 9, 10, 11, 11, 12, 13, 13, 14, 13,17];
pixel_mat_y = [19, 17, 16, 15, 14, 13, 12, 11, 10, 9, 10, 8, 6, 5, 4, 5, 4, 4, 3, 3, 3, 3, 3, 4, 4, 6, 9, 9, 10, 11, 12, 14, 15, 16, 17, 15, 17, 18, 18, 20, 18, 18, 18, 19, 20, 20,19];
num_of_pixels = length(pixel_mat_y);

% number of the angle bins
num_of_angle_bins = 58;
% num_of_angle_bins = 360;

% Define the center of the ring and outer radius
x_center = 12;
y_center = 12;
outer_radius = 9.5;
outer_radius_error = 0.025;

% matrix to save the correlation graphs versus angles
correlation_vecs_linear = zeros(num_of_pixels, num_of_angle_bins);

% for n=1:1
for n=1:num_of_pixels
    
    x_coord = pixel_mat_x(n);
    y_coord = pixel_mat_y(n);
    
    y_original=reshape(cor_matS_C_linear(x_coord,y_coord,:,:), 26, 26);
    y = y_original;
    
    % zero out pixels near the origin
    for i=-4:1:4
        for j=-4:1:4
  
            if x_coord+i > 26 || y_coord+j > 26 || x_coord+i < 1 || y_coord+j < 1
                continue
            end

            y(x_coord+i, y_coord+j) = 1;

        end
    end

    z=real(log(y));
    t=sign(y);
            
    % zero out mask with negative entries
    for i=1:26

        for j=1:26
                    
            if sign(t(i, j)) < 0
                t(i, j) = 0;
            end
        end
    end


    % this is current correlation image
    current_corr_image = z.*t.*fil;

    % get the point corresponding to the reference pixel on the outer
    % radius 
    out_ring_coords = find_point_on_outer_radius([x_center, y_center], [x_coord y_coord], outer_radius);

    % now get the set of points that goes over the outer ring starting from
    % the reference pixel
    out_ring_list = generate_points_on_outer_radius_clockwise([x_center, y_center], [out_ring_coords(1), out_ring_coords(2)], outer_radius, num_of_angle_bins);


    % we search over the whole ring given starting point above which is (x_coord,
    % y_coord)
    for ring_point_id=1:length(out_ring_list)

        x_ring_point = out_ring_list(ring_point_id, 1);
        y_ring_point = out_ring_list(ring_point_id, 2);


        % check if we are indeed on the outer ring, since the grid is very
        % sparse, we allow for some small error with the location on the
        % outer ring
        if (x_ring_point - x_center)^2 + (y_ring_point - y_center)^2 >= (outer_radius)^2 || ...
           (x_ring_point - x_center)^2 + (y_ring_point - y_center)^2 <= (outer_radius + outer_radius_error)^2
                
            % get all the points on the line between the outer ring point
            % and the center
            points_from_center_to_outer_ring = bresenham_line([x_center, y_center], [x_ring_point, y_ring_point]);
            
            
            % go over the line and check if we have some non-zero value
            for line_point_id=1:length(points_from_center_to_outer_ring)
    
                x_line_point = points_from_center_to_outer_ring(line_point_id, 1);
                y_line_point = points_from_center_to_outer_ring(line_point_id, 2);
                
                if current_corr_image(x_line_point, y_line_point) > 0

                    correlation_vecs_linear(n, ring_point_id) = 1;
                    (ring_point_id / 58) * 360
                    break

    
                end 
            end
        end            
    end 

end

% matrix to save the correlation graphs versus angles
correlation_vecs_circular = zeros(num_of_pixels, num_of_angle_bins);

% for n=1:1
for n=1:num_of_pixels
    
    x_coord = pixel_mat_x(n);
    y_coord = pixel_mat_y(n);
    
    y_original=reshape(cor_matS_C_circular(x_coord,y_coord,:,:), 26, 26);
    y = y_original;
    
    % zero out pixels near the origin
    for i=-4:1:4
        for j=-4:1:4

            if x_coord+i > 26 || y_coord+j > 26 || x_coord+i < 1 || y_coord+j < 1
                continue
            end

            y(x_coord+i, y_coord+j) = 1;

        end
    end

    z=real(log(y));
    t=sign(y);
            
    % zero out mask with negative entries
    for i=1:26

        for j=1:26
                    
            if sign(t(i, j)) < 0
                t(i, j) = 0;
            end
        end
    end


    % this is current correlation image
    current_corr_image = z.*t.*fil;

    % get the point corresponding to the reference pixel on the outer
    % radius 
    out_ring_coords = find_point_on_outer_radius([x_center, y_center], [x_coord y_coord], outer_radius);

    % now get the set of points that goes over the outer ring starting from
    % the reference pixel
    out_ring_list = generate_points_on_outer_radius_clockwise([x_center, y_center], [out_ring_coords(1), out_ring_coords(2)], outer_radius, num_of_angle_bins);


    % we search over the whole ring given starting point above which is (x_coord,
    % y_coord)
    for ring_point_id=1:length(out_ring_list)

        x_ring_point = out_ring_list(ring_point_id, 1);
        y_ring_point = out_ring_list(ring_point_id, 2);

        % check if we are indeed on the outer ring, since the grid is very
        % sparse, we allow for some small error with the location on the
        % outer ring
        if (x_ring_point - x_center)^2 + (y_ring_point - y_center)^2 >= (outer_radius)^2 || ...
           (x_ring_point - x_center)^2 + (y_ring_point - y_center)^2 <= (outer_radius + outer_radius_error)^2
                
            % get all the points on the line between the outer ring point
            % and the center
            points_from_center_to_outer_ring = bresenham_line([x_center, y_center], [x_ring_point, y_ring_point]);
            
            
            % go over the line and check if we have some non-zero value
            for line_point_id=1:length(points_from_center_to_outer_ring)
    
                x_line_point = points_from_center_to_outer_ring(line_point_id, 1);
                y_line_point = points_from_center_to_outer_ring(line_point_id, 2);
                
                if current_corr_image(x_line_point, y_line_point) > 0

                    correlation_vecs_circular(n, ring_point_id) = 1;
                    (ring_point_id / 58) * 360
                    break

    
                end 
            end
        end            
    end 

end


corr_sum_circular = sum(correlation_vecs_circular, 1);
corr_sum_linear = sum(correlation_vecs_linear, 1);
% corr_mat_sum = zeros(58, 2);
% for i=1:1:58
% 
%     corr_mat_sum(i, 1) = corr_sum_circular(i) / sum(corr_sum_circular);
%     corr_mat_sum(i, 2) = corr_sum_linear(i) / sum(corr_sum_linear);
% 
% end


bin_vec = [1:1:58];


figure()
% bar(bin_vec, corr_mat_sum)
bar_color = [146 18 39] / 256; % Blue color
alpha_value = 0.6;
bar_circular = bar(bin_vec, corr_sum_circular / sum(corr_sum_circular))
% Set the face color and transparency for each bar
for i = 1:numel(bar_circular)
    set(bar_circular(i), 'FaceColor', bar_color, 'FaceAlpha', alpha_value);
end
hold on


% bar(bin_vec, corr_mat_sum)
bar_color = [18 65 146] / 256; % greenish color
alpha_value = 0.6;
bar_linear = bar(bin_vec, corr_sum_linear / sum(corr_sum_linear))
% Set the face color and transparency for each bar
for i = 1:numel(bar_linear)
    set(bar_linear(i), 'FaceColor', bar_color, 'FaceAlpha', alpha_value);
end
axis tight


curve_to_fit = corr_sum_linear / sum(corr_sum_linear);


function points_on_line = bresenham_line(center_coords, outer_point_coords)
    % Extract coordinates of points V and W
    x1 = center_coords(1);
    y1 = center_coords(2);
    x2 = outer_point_coords(1);
    y2 = outer_point_coords(2);
    
    % Calculate differences and absolute differences
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    
    % Determine the sign of increments
    if x1 < x2
        sx = 1;
    else
        sx = -1;
    end
    
    if y1 < y2
        sy = 1;
    else
        sy = -1;
    end
    
    % Initialize error and error increments
    err = dx - dy;
    
    % Initialize coordinates
    x = x1;
    y = y1;
    
    % Initialize list to store points on the line
    points_on_line = [];
    
    % Iterate over the line using Bresenham's algorithm
    while true
        % Add current point to the list
        points_on_line = [points_on_line; [x, y]];
        
        % Check if the end point is reached
        if x == x2 && y == y2

            break;
        end
        
        % Calculate error for next step
        e2 = 2 * err;
        if e2 > -dy
            err = err - dy;
            x = x + sx;
        end
        
        if e2 < dx
            err = err + dx;
            y = y + sy;
        end
    end
end




function point_on_outer_radius = find_point_on_outer_radius(center_coords, reference_point_coords, outer_radius)
    % Calculate polar coordinates of point V relative to the center
    [theta, ~] = cart2pol(reference_point_coords(1) - center_coords(1), reference_point_coords(2) - center_coords(2));
    
    % Calculate the angle between V and the positive x-axis
%     theta_outer = rad2deg(theta);
    
    % Calculate the angle between the positive x-axis and the line connecting
    % the center of the ring and the point on the outer radius
%     angle_outer = angle + 180; % Add 180 degrees to obtain the angle on the outer radius
    
    % Convert the angle to radians
%     theta_outer = deg2rad(angle_outer);
    
    % Calculate the coordinates of the point on the outer radius with
    % respect to (0,0) position!!!
    x_outer = center_coords(2) + outer_radius * sin(theta);
    y_outer = center_coords(1) + outer_radius * cos(theta);

    % let is round up the values according to the quadrant we are in
    x_outer_rounded = 0;
    y_outer_rounded = 0;

    if theta > 0 && theta <= 90

        x_outer_rounded = ceil(x_outer);
        y_outer_rounded = floor(y_outer);

    elseif theta > 90 && theta <= 180

        x_outer_rounded = ceil(x_outer);
        y_outer_rounded = ceil(y_outer);

    elseif theta > 180 && theta <= 270

        x_outer_rounded = floor(x_outer);
        y_outer_rounded = ceil(y_outer);

    else

        x_outer_rounded = floor(x_outer);
        y_outer_rounded = floor(y_outer);

    end 



    
    % Return the coordinates of the point on the outer radius and don't
    % forget to round them up to the closest integer from above
    point_on_outer_radius = [x_outer_rounded, y_outer_rounded];
end


function points_on_outer_radius_clockwise = generate_points_on_outer_radius_clockwise(center_coords, outer_ring_point, outer_radius, num_of_angle_bins)
    % Calculate polar coordinates of point P relative to the center
    [theta_P, ~] = cart2pol(outer_ring_point(1) - center_coords(1), outer_ring_point(2) - center_coords(2));
    
    % Initialize list to store points on the outer radius
    points_on_outer_radius_clockwise = [];
    
    % Iterate over angles from 0 to 2*pi radians
    for theta = linspace(theta_P + (pi / 6), theta_P + 2*pi - (pi / 6), num_of_angle_bins)
        % Calculate the coordinates of the point on the outer radius
        x = center_coords(2) + outer_radius * sin(theta);
        y = center_coords(1) + outer_radius * cos(theta);

        % let is round up the values according to the quadrant we are in
        x_rounded = 0;
        y_rounded = 0;

        if theta > 0 && theta <= 90

            x_rounded = ceil(x);
            y_rounded = floor(y);

        elseif theta > 90 && theta <= 180

            x_rounded = ceil(x);
            y_rounded = ceil(y);

        elseif theta > 180 && theta <= 270

            x_rounded = floor(x);
            y_rounded = ceil(y);

        else

            x_rounded = floor(x);
            y_rounded = floor(y);

        end 
        
        % Add the point to the list and don't forget to round them to the
        % closest integer from below
        points_on_outer_radius_clockwise = [points_on_outer_radius_clockwise; [x_rounded, y_rounded]];
    end
end





