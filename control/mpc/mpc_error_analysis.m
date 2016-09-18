% mpc_error_analysis.m
% Copyright 2016 Andrew P. Sabelhaus, Berkeley Emergent Space Tensegrities Lab
% This function loads in a set of saved data from an individual MPC run, calculates
% the error in trajectory tracking, and ...

function [ errors ] = mpc_error_analysis( start_time_string, path_to_data_folder, plots_flag )
% Inputs:
%   start_time_string = the prefix of the data file to read in, for analysis.
%       note that the remainder of the string is hard-coded below. TO-DO: change this.
%   path_to_data_folder = location of the data file to read in.
%       NOTE that this script currently ONLY WORKS for three-vertebrae spines.
%   plots_flag = create plots (1) or do not make plots (0).
% Outputs:
%   errors = a struct containing all the types of errors that were calculated 

% Some hard-coded variables 
% TO-DO: pass these in to the function?
% Size of the text in the figures
fontsize = 12;

% Read in the data file.
data_path = strcat( path_to_data_folder, 'ultra-spine-mpc_data_', start_time_string );
data = load(data_path);

% A quick sanity check: only 3-vertebra trajectories allowed as of 2016-04-25.
assert( size(data.traj, 1) == 36, 'MPC error script only works for three-vertebra spines. Use a different data file.');

% Pick out the reference trajectory and the resulting trajectory.
% Each of these should be num_states x length_trajectory
ref_traj = data.traj;
result_traj = data.refx;

% ref_traj is off by one. Prune out the first point in the reference trajectory (this point, the intial state, is not controlled for.)
ref_traj = ref_traj(:, 2:end);

% Append these to the errors struct.
errors.ref_traj = ref_traj;
errors.result_traj = result_traj;

% Take the difference between these two trajectories, and square each point: this is the squared tracking error (all states)
tracking_error_squared = (result_traj - ref_traj).^2;
% what does just the value of the error look like? These units will be in meters, at least, for comparison.
tracking_error = result_traj - ref_traj;
% Save this into the output struct
errors.tracking_error_squared = tracking_error_squared;
errors.tracking_error = tracking_error;

% Take the total tracking error, as a function of timestep, for the trajectories that were actually tracked
% using MPC. (ie., remove the tracking for velocities.)
% tracking_XYZ will be for the positions of each tetrahedron. 3 tetras, 3 states, = 9 trajectories.
% tracking_angle will be for rotations
tracking_XYZ_squared = zeros(9, size(tracking_error_squared, 2));
tracking_angle_squared = zeros(9, size(tracking_error_squared, 2));
tracking_XYZ = zeros(9, size(tracking_error, 2));
tracking_angle = zeros(9, size(tracking_error, 2));

for i=1:3
    % Pick out the XYZ and angle states
    tracking_XYZ_squared( (i-1)*3 + 1: (i-1)*3 + 3, :) = tracking_error_squared( (i-1)*12 + 1: (i-1)*12 + 3, :);
    tracking_angle_squared( (i-1)*3 + 1: (i-1)*3 + 3, :) = tracking_error_squared( (i-1)*12 + 4: (i-1)*12 + 6, :);
    tracking_XYZ( (i-1)*3 + 1: (i-1)*3 + 3, :) = tracking_error( (i-1)*12 + 1: (i-1)*12 + 3, :);
    tracking_angle( (i-1)*3 + 1: (i-1)*3 + 3, :) = tracking_error( (i-1)*12 + 4: (i-1)*12 + 6, :);
end

% save these into the output struct
errors.tracking_XYZ_squared = tracking_XYZ_squared;
errors.tracking_angle_squared = tracking_angle_squared;
errors.tracking_XYZ = tracking_XYZ;
errors.tracking_angle = tracking_angle;

% Smash each of these down to a single time series of scalars
tracking_XYZ_squared_sum = sum(tracking_XYZ_squared);
tracking_angle_squared_sum = sum(tracking_angle_squared);
tracking_XYZ_sum = sum(tracking_XYZ);
tracking_angle_sum = sum(tracking_angle);

% save these into the output struct
errors.tracking_XYZ_squared_sum = tracking_XYZ_squared_sum;
errors.tracking_angle_squared_sum = tracking_angle_squared_sum;
errors.tracking_XYZ_sum = tracking_XYZ_sum;
errors.tracking_angle_sum = tracking_angle_sum;

% add these two to get the total error for all position states
tracking_squared_total = tracking_XYZ_squared_sum + tracking_angle_squared_sum;
errors.tracking_squared_total = tracking_squared_total;
tracking_total = tracking_XYZ_sum + tracking_angle_sum;
errors.tracking_total = tracking_total;

% Make plots, if flag is set
if plots_flag
    % Calculate the time in seconds at each point.
    t = data.dt: data.dt : (data.num_points_ref_traj-1)*data.dt;
    % Adjust time so that everything is in milliseconds. This makes things clearer.
    t = t * 100;
    % A good figure window setup is 'Position',[100,100,500,300].
    size(t)
    
    % Make a sequence of points at zero for plotting a thin dashed black line in the below plots.
    %zero_line = zeros(1, size(t,2));
    %size(zero_line)
    
    %% Plot the errors, all of XYZTGP on the same subplot.
    % Create the handle for the overall figure
    % Also, use the OpenGL renderer so that symbols are formatted correctly.
    %errors_handle = figure('Renderer', 'opengl');
    errors_handle = figure;
    hold on;
    set(gca,'FontSize',fontsize);
    % This figure will have 6 smaller plots, so make it twice the size of my usual window dimensions.
    % Actually, since this will be 3x2, needs to be 2x in the lateral direction.
    %set(errors_handle,'Position',[100,100,450,600]);
    set(errors_handle,'Position',[100,100,900,300]);
    % Start the first subplot
    subplot_handle = subplot(3, 2, 1);
    hold on;
    
    % Plot the X errors
    for i=1:3
        % Scale the lengths here to get centimeters.
        plot(t, tracking_XYZ( (i-1) * 3 + 1, :)*100, 'LineWidth',2 );
    end
    % Plot the zero line
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %xlabel('Time (msec)');
    ylabel('e_X (cm)');
    % Only create a title for the first plot, that will serve for all the others too.
    %title('Tracking Errors in X Y Z  \theta \gamma \psi');
    title('   Postitions (No Disturbance)');
    % Scale the plot. A good scale here is...
    %ylim([-1.0 1.5]);
    ylim([-2.0 3.5]);
    
    % Move the plot very slightly to the left
    P = get(subplot_handle,'Position')
    %set(subplot_handle,'Position',[P(1)-0.05 P(2) P(3) P(4)])
    set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.05 P(3)+0.01 P(4)-0.04])
    hold off;
    
    % Title the whole plot
    %mtit('Tracking Errors, With Disturbances','xoff',0.5,'yoff',0.2)
    

    % Plot the Y errors
    subplot_handle = subplot(3, 2, 3);
    hold on;
    for i=1:3
        % Scale the lengths here to get centimeters.
        plot(t, tracking_XYZ( (i-1) * 3 + 2, :)*100, 'LineWidth',2  );
    end
    % Plot the zero line
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    %xlabel('Time (msec)');
    ylabel('e_Y (cm)');
    %title('Tracking Error in Y');
    % Scale the plot. A good scale here is...
    % used to be -0.2 to 0.1
    %ylim([-1 1]);
    ylim([-2 3.5]);
    %ylim([-2 2]);    
    % Move the plot very slightly to the left
    % For these lower figures, move them upwards a bit more.
    P = get(subplot_handle,'Position')
    set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.07 P(3)+0.01 P(4)-0.04])
    
    hold off;
    
    
    % Plot the Z errors
    subplot_handle = subplot(3, 2, 5);
    hold on;
    for i=1:3
        % Scale the lengths here to get centimeters.
        plot(t, tracking_XYZ( (i-1) * 3 + 3, :)*100, 'LineWidth',2 );
    end
    % Plot the zero line
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    %xlabel('Time (msec)');
    ylabel('e_Z (cm)');
    %title('Tracking Error in Z');
    % Scale the plot. A good scale here is...
    % Used to be -0.2 to 0.2
    %ylim([-1 1]);
    ylim([-2 3.5]);
    %ylim([-2 2]);
    % make a x-label here, for the position states.
    xlabel('Time (msec)');    
    % Move the plot very slightly to the left
    P = get(subplot_handle,'Position')
    set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.10 P(3)+0.01 P(4)-0.04])
    
    hold off;
    
    % Plot the Theta errors:
    subplot_handle = subplot(3, 2, 2);
    hold on;
    for i=1:3
        % Convert these radians to degrees.
        plot(t, tracking_angle( (i-1) * 3 + 1, :) * (180/pi), 'LineWidth',2 );
    end
    % Plot the zero line
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    %xlabel('Time (msec)');
    %ylabel('\theta (deg)');
    %ylabel('e_T (deg,\circ)');
    ylabel('e_\theta (deg)');
    %title('Tracking Error in \theta');
    % Scale the plot. A good scale here is...
    % Used to be: [-1 1]
    ylim([-10 15]);
    %ylim([-2 2]);
    % Make a title here, for the angle states
    title('   Angles (No Disturbance)');    
    % Move the plot very slightly to the left
    P = get(subplot_handle,'Position')
    set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.05 P(3)+0.01 P(4)-0.04])
    
    hold off;
    
    % Plot the Gamma errors:
    subplot_handle = subplot(3, 2, 4);
    hold on;
    for i=1:3
        % Convert these radians to degrees.
        plot(t, tracking_angle( (i-1) * 3 + 2, :) * (180/pi), 'LineWidth',2 );
    end
    % Plot the zero line
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    %xlabel('Time (msec)');
    %ylabel('\gamma (deg)');
    ylabel('e_\gamma (deg)');
    %title('Tracking Error in \gamma');
    % Scale the plot. A good scale here is...
    ylim([-10 15]);    
    % Move the plot very slightly to the left
    P = get(subplot_handle,'Position')
    set(subplot_handle,'Position',[[P(1)-0.06 P(2)+0.07 P(3)+0.01 P(4)-0.04]])
    
    hold off;
    
    % Plot the Psi errors:
    subplot_handle = subplot(3, 2, 6);
    hold on;
    for i=1:3
        % Convert these radians to degrees.
        plot(t, tracking_angle( (i-1) * 3 + 3, :) * (180/pi), 'LineWidth',2 );
    end
    % Plot the zero line
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    xlabel('Time (msec)');
    %ylabel('\psi (deg)');
    ylabel('e_\psi (deg)');    
    % Move the plot very slightly to the left
    P = get(subplot_handle,'Position')
    set(subplot_handle,'Position',[[P(1)-0.06 P(2)+0.10 P(3)+0.01 P(4)-0.04]])
    
    
    % Make the legend
    % Item names: use newlines.
    vert1label = sprintf('Vertebra 1 \n (Bottom)');
    vert2label = sprintf('Vertebra 2 \n (Middle)');
    vert3label = sprintf('Vertebra 3 \n (Top)');
    legend_handle = legend(vert1label, vert2label, vert3label);
    % move the legend to the corner of the figure.
    % Position is: [left, bottom, width, height]
    % These units are expressed as percents of the total figure
    % Bottom left, for the 1x6 plot:
    %legend_position = [0.22 0.042 0 0];
    % Bottom right, 1x6 plot:
    %legend_position = [0.82 0.042 0 0];
    % Bottom right, 3x2 plot:
    legend_position = [0.935 0.6 0.001 0];
    set(legend_handle,'Position', legend_position);
    %title('Tracking Error in \psi');
    % Scale the plot. A good scale here is...
    % Used to be: [-2 2]
    ylim([-10 15]);
    %ylim([-2 2]);
    hold off;
    
    % some debugging to figure out how to shift the subplot down a bit and add an overall title
    %drawnow;
    %disp('subplot handle position:');
    %P = get(subplot_handle,'Position')
    %set(subplot_handle,'Position',[P(1)+0.1 P(2) P(3) P(4)]);
    % ...didn't work, just did separate labels.
    
    %% Plot the XZ positions of the reference trajectory and resulting trajectory.
    % This should give a good visual representation of the spine's movement.
    
    % Top:
    trajectories_top_handle = figure;
    hold on;
    set(trajectories_top_handle,'Position',[100,100,500,300]);
    % Plot the reference: X vs. Z. The top vertebra is at states 25 to 36.
    % Scale the lengths here to get cm.
    plot( ref_traj(25,:)*100, ref_traj(27, :)*100, 'b.');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj(25,:)*100, result_traj(27,:)*100, 'g.');
    legend('Reference', 'Result');
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of Top Vertebra');
    % Scale the plot?
    hold off;
    
    % Middle:
    trajectories_middle_handle = figure;
    hold on;
    set(trajectories_middle_handle,'Position',[100,100,500,300]);
    % Plot the reference: X vs. Z. The middle vertebra is at states 13 to 24.
    % Scale the lengths here to get cm.
    plot( ref_traj(13,:)*100, ref_traj(15, :)*100, 'b.');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj(13,:)*100, result_traj(15,:)*100, 'g.');
    legend('Reference', 'Result');
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of Middle Vertebra');
    % Scale the plot?
    hold off;
    
    % Lower:
    trajectories_lower_handle = figure;
    hold on;
    set(trajectories_lower_handle,'Position',[100,100,500,300]);
    % Plot the reference: X vs. Z. The lower vertebra is at states 1 to 12.
    % Scale the lengths here to get cm.
    plot( ref_traj(1,:)*100, ref_traj(3, :)*100, 'b.');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj(1,:)*100, result_traj(3,:)*100, 'g.');
    legend('Reference', 'Result');
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of Lower Vertebra');
    % Scale the plot?
    hold off;
    
    %% Plot the combined position of all three vertebrae on the same plot.
    % Will this be interesting??
    
    % Create the figure handle: make it large
    % Also, use the OpenGL renderer so that symbols are formatted correctly.
    all_positions_handle = figure('Renderer', 'opengl');
    hold on;
    set(gca,'FontSize',fontsize);
    % This figure will have 6 smaller plots, so make it twice the size of my usual window dimensions.
    set(all_positions_handle,'Position',[100,100,500,300]);
    % Set some reasonable limits.
    ylim([9 31]);
    xlim([-12 0]);
    
    % Title, label, etc
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of vertebrae during tracking control');
    
    % Plot the reference: X vs. Z. The top vertebra is at states 25 to 36.
    % Scale the lengths here to get cm.
    plot( ref_traj(25,:)*100, ref_traj(27, :)*100, 'b.');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj(25,:)*100, result_traj(27,:)*100, 'g.');
    
    % Plot the reference: X vs. Z. The middle vertebra is at states 13 to 24.
    % Scale the lengths here to get cm.
    plot( ref_traj(13,:)*100, ref_traj(15, :)*100, 'b.');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj(13,:)*100, result_traj(15,:)*100, 'g.');
    
    % Plot the reference: X vs. Z. The lower vertebra is at states 1 to 12.
    % Scale the lengths here to get cm.
    plot( ref_traj(1,:)*100, ref_traj(3, :)*100, 'b.');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj(1,:)*100, result_traj(3,:)*100, 'g.');
    
    % Make a legend. Since all vertebrae are the same color, we only need two labels.
    legend('Reference Trajectories, in X-Z', 'Result Traj., No Noise', 'Location', 'Southwest');
    
    %% Plot the total sum-squared error for positions and angles
    
    % Since it doesn't make sense to combine errors, since the magnitudes are different,
    % provide 2 different plots.
    % Also, use the OpenGL renderer so that symbols are formatted correctly.
    total_error_handle = figure('Renderer', 'opengl');
    hold on;
    set(total_error_handle,'Position',[100,100,500,300]);
    % Make a subplot:
    subplot(2, 1, 1);
    % Adjust the error from meters to centimeters: since it's squared,
    % Changing m^2 to cm^2 requires mulitplication by 100^2
    plot( t, tracking_XYZ_squared_sum, '.-');
    xlabel('Time (msec)');
    ylabel( sprintf('Sq. Err. (m^2)') );
    title('Total squared error for XYZ states')
    % Scale this plot to emphasize how small these errors are.
    ylim([0 0.00125]);
    % Make the font larger for these subplots that get squished.
    %set(gca,'FontSize',13);
    % Make the second plot:
    subplot(2, 1, 2);
    % Adjust this error to degrees squared, since that's more intuitive than rad^2.
    % That means multiply by (180/pi)^2
    plot(t, tracking_angle_squared_sum, '.-');
    xlabel('Time (msec)');
    ylabel( sprintf('Sq. Err. (rad^2)') );
    title('Total squared error for angle (\theta \gamma \psi) states');
    % Make the font larger for these subplots that get squished.
    %set(gca,'FontSize',13);
    % Scale:
    ylim([0 0.05]);
    hold off;
    
    % Let's try to combine anyway. Normalize by dividing by the mean value in each vector.
    %total_combined_error_handle = figure;
    %hold on;
    %set(total_combined_error_handle,'Position',[100,100,500,300]);
    % The mean value of the XYZ states:
    %m_XYZ = mean(tracking_XYZ_squared_sum);
    % The mean value of the TGP states:
    %m_TGP = mean(tracking_angle_squared_sum);
    % calculate the normalized errors:
    %XYZ_norm = tracking_XYZ_squared_sum - m_XYZ;
    %TGP_norm = tracking_angle_squared_sum - m_TGP;
    % Add them and plot:
    %total_combined_error = XYZ_norm + TGP_norm;
    %hold off;
    %plot(t, XYZ_norm);
    %title('XYZ_norm');
    %figure;
    %plot(t, TGP_norm);
    %title('TGP_norm');
    %figure;
    %plot(t, total_combined_error);
    %title('total_combined_error');
    
    %% Calculate the "max errors" as a metric to use:
    %max_error_x = 

end

end






