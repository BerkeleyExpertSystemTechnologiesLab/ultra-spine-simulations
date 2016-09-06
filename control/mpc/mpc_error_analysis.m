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
    % Plot the errors in X,Y,Z: absolute error for each
    
    %% Plot the errors. First, XYZ positions:
    % Plot the X errors
    errors_X_figure_handle = figure;
    hold on;
    set(gca,'FontSize',fontsize);
    set(errors_X_figure_handle,'Position',[100,100,500,300]);
    for i=1:3
        % Scale the lengths here to get centimeters.
        plot(t, tracking_XYZ( (i-1) * 3 + 1, :)*100 );
    end
    legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    xlabel('Time (msec)');
    ylabel('Error in X (cm)');
    title('Tracking Error in X');
    % Scale the plot. A good scale here is...
    ylim([-2 2]);
    hold off;
    
    % Plot the Y errors
    errors_Y_figure_handle = figure;
    hold on;
    set(gca,'FontSize',fontsize);
    set(errors_Y_figure_handle,'Position',[100,100,500,300]);
    for i=1:3
        % Scale the lengths here to get centimeters.
        plot(t, tracking_XYZ( (i-1) * 3 + 2, :)*100 );
    end
    legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    xlabel('Time (msec)');
    ylabel('Error in Y (cm)');
    title('Tracking Error in Y');
    % Scale the plot. A good scale here is...
    ylim([-2 2]);
    hold off;
    
    % Plot the Z errors
    errors_Z_figure_handle = figure;
    hold on;
    set(gca,'FontSize',fontsize);
    set(errors_Z_figure_handle,'Position',[100,100,500,300]);
    for i=1:3
        % Scale the lengths here to get centimeters.
        plot(t, tracking_XYZ( (i-1) * 3 + 3, :)*100 );
    end
    legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    xlabel('Time (msec)');
    ylabel('Error in Z (cm)');
    title('Tracking Error in Z');
    % Scale the plot. A good scale here is...
    ylim([-2 2]);
    hold off;
    
    %% Plot the errors in TGP:
    
    % Plot the Theta errors:
    errors_T_figure_handle = figure;
    hold on;
    set(gca,'FontSize',fontsize);
    set(errors_T_figure_handle,'Position',[100,100,500,300]);
    for i=1:3
        % Convert these radians to degrees.
        plot(t, tracking_angle( (i-1) * 3 + 1, :) * (180/pi) );
    end
    legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    xlabel('Time (msec)');
    ylabel('Error in \theta (degrees)');
    title('Tracking Error in \theta');
    % Scale the plot. A good scale here is...
    %ylim([-2 2]);
    hold off;
    
    % Plot the Gamma errors:
    errors_G_figure_handle = figure;
    hold on;
    set(gca,'FontSize',fontsize);
    set(errors_G_figure_handle,'Position',[100,100,500,300]);
    for i=1:3
        % Convert these radians to degrees.
        plot(t, tracking_angle( (i-1) * 3 + 2, :) * (180/pi) );
    end
    legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    xlabel('Time (msec)');
    ylabel('Error in \gamma (degrees)');
    title('Tracking Error in \gamma');
    % Scale the plot. A good scale here is...
    %ylim([-2 2]);
    hold off;
    
    % Plot the Psi errors:
    errors_P_figure_handle = figure;
    hold on;
    set(gca,'FontSize',fontsize);
    set(errors_P_figure_handle,'Position',[100,100,500,300]);
    for i=1:3
        % Convert these radians to degrees.
        plot(t, tracking_angle( (i-1) * 3 + 3, :) * (180/pi) );
    end
    legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    xlabel('Time (msec)');
    ylabel('Error in \psi (degrees)');
    title('Tracking Error in \psi');
    % Scale the plot. A good scale here is...
    %ylim([-2 2]);
    hold off;
    
    
    
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
    ylim([0 0.001]);
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
    
    % Plot the error in position for all three of the vertebrae together

end

end






