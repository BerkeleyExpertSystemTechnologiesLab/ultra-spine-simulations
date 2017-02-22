% mpc_error_analysis_combined.m
% Copyright 2016 Andrew P. Sabelhaus, Berkeley Emergent Space Tensegrities Lab
% This function loads in a set of saved data from multiple MPC runs and plots a comparison

function mpc_error_analysis_combined( start_time_string1, start_time_string2, path_to_data_folder, plots_flag )
% Inputs:
%   start_time_string1 = the prefix of the data file to read in, for analysis.
%       note that the remainder of the string is hard-coded below. TO-DO: change this.
%   start_time_string2 = prefix of the second set of data
%   path_to_data_folder = location of the data file to read in.
%       NOTE that this script currently ONLY WORKS for three-vertebrae spines.
%   plots_flag = create plots (1) or do not make plots (0).
% Outputs:
%   none, since this function just plots.

% Some hard-coded variables 
% TO-DO: pass these in to the function?
% Size of the text in the figures
fontsize = 12;

% Call mpc_error_analysis to calculate all the errors for both MPC runs.
errors1 = mpc_error_analysis( start_time_string1, path_to_data_folder, 0);
errors2 = mpc_error_analysis( start_time_string2, path_to_data_folder, 0);

% Make plots, if flag is set
if plots_flag
    % Calculate the time in seconds at each point.
    % Assume that the times for errors1 are the same as for errors2.
    t = errors1.dt: errors1.dt : (errors1.num_points_ref_traj-1)*errors1.dt;
    % Adjust time so that everything is in milliseconds. This makes things clearer.
    t = t * 100;
    % A good figure window setup is 'Position',[100,100,500,300].
    size(t)
    
    % Pick out some of the variables so that it's easier to write the code below:
    % Assume that both the data files use the same reference trajectory.
    ref_traj = errors1.ref_traj;
    result_traj1 = errors1.result_traj;
    result_traj2 = errors2.result_traj;
    tracking_XYZ_squared_sum1 = errors1.tracking_XYZ_squared_sum;
    tracking_XYZ_squared_sum2 = errors2.tracking_XYZ_squared_sum;
    tracking_angle_squared_sum1 = errors1.tracking_angle_squared_sum;
    tracking_angle_squared_sum2 = errors2.tracking_angle_squared_sum;
    
    %% Plot the XZ positions of the reference trajectory and resulting trajectory.
    % This should give a good visual representation of the spine's movement.
    
    % Top:
    trajectories_top_handle = figure;
    hold on;
    set(trajectories_top_handle,'Position',[100,100,500,300]);
    % Set the axis limits
    xlim([-12.5 0.1])
    ylim([27.6 30.1]);
    % Plot the reference: X vs. Z. The top vertebra is at states 25 to 36.
    % Scale the lengths here to get cm.
    plot( ref_traj(25,:)*100, ref_traj(27, :)*100, 'b','LineWidth',2);
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj1(25,:)*100, result_traj1(27,:)*100, 'g','LineWidth',2);
    plot( result_traj2(25,:)*100, result_traj2(27,:)*100, 'm','LineWidth',2);
    legend('Reference', 'Result, No Dist.','Result, With Dist.','Location','Southeast');
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of Top Vertebra');
    set(gca,'FontSize',fontsize);
    % Scale the plot?
    hold off;
    
    % Middle:
    trajectories_middle_handle = figure;
    hold on;
    set(trajectories_middle_handle,'Position',[100,100,500,300]);
    % Plot the reference: X vs. Z. The middle vertebra is at states 13 to 24.
    % Scale the lengths here to get cm.
    plot( ref_traj(13,:)*100, ref_traj(15, :)*100, 'b','LineWidth',2);
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj1(13,:)*100, result_traj1(15,:)*100, 'g','LineWidth',2);
    plot( result_traj2(13,:)*100, result_traj2(15,:)*100, 'm','LineWidth',2);
    legend('Reference', 'Result, No Dist.', 'Result, With Dist.','Location','Southeast');
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of Middle Vertebra');
    set(gca,'FontSize',fontsize);
    % Scale the plot?
    hold off;
    
    % Lower:
    trajectories_lower_handle = figure;
    hold on;
    set(trajectories_lower_handle,'Position',[100,100,500,300]);
    % Plot the reference: X vs. Z. The lower vertebra is at states 1 to 12.
    % Scale the lengths here to get cm.
    plot( ref_traj(1,:)*100, ref_traj(3, :)*100, 'b','LineWidth',2);
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj1(1,:)*100, result_traj1(3,:)*100, 'g','LineWidth',2);
    plot( result_traj2(1,:)*100, result_traj2(3,:)*100, 'm','LineWidth',2);
    legend('Reference', 'Result, No Dist.', 'Result, With Dist.','Location','Southeast');
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of Lower Vertebra');
    set(gca,'FontSize',fontsize);
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
    plot( ref_traj(25,:)*100, ref_traj(27, :)*100, 'b.-');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj1(25,:)*100, result_traj1(27,:)*100, 'g.-');
    plot( result_traj2(25,:)*100, result_traj1(27,:)*100, 'm.-');
    
    % Plot the reference: X vs. Z. The middle vertebra is at states 13 to 24.
    % Scale the lengths here to get cm.
    plot( ref_traj(13,:)*100, ref_traj(15, :)*100, 'b.-');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj1(13,:)*100, result_traj1(15,:)*100, 'g.-');
    plot( result_traj2(13,:)*100, result_traj1(15,:)*100, 'm.-');
    
    % Plot the reference: X vs. Z. The lower vertebra is at states 1 to 12.
    % Scale the lengths here to get cm.
    plot( ref_traj(1,:)*100, ref_traj(3, :)*100, 'b.-');
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    plot( result_traj1(1,:)*100, result_traj1(3,:)*100, 'g.-');
    plot( result_traj2(1,:)*100, result_traj1(3,:)*100, 'm.-');
    
    % Make a legend. Since all vertebrae are the same color, we only need two labels.
    legend('Reference Trajs.', 'Result, No Dist.', 'Result, With Dist.', 'Location', 'Southwest');
    
    %% Plot the total sum-squared error for positions and angles
    
    % Since it doesn't make sense to combine errors, since the magnitudes are different,
    % provide 2 different plots.
    % Also, use the OpenGL renderer so that symbols are formatted correctly.
    %total_error_handle = figure('Renderer', 'opengl');
    total_error_handle = figure;
    hold on;
    set(total_error_handle,'Position',[100,100,500,300]);
    % Make a subplot:
    subplot(2, 1, 1);
    hold on;
    % Adjust the error from meters to centimeters: since it's squared,
    % Changing m^2 to cm^2 requires mulitplication by 100^2
    plot( t, sqrt(tracking_XYZ_squared_sum1)*100, 'g.-');
    plot( t, sqrt(tracking_XYZ_squared_sum2)*100, 'm.-');
    xlabel('Time (msec)');
    ylabel( sprintf('Sum Err. (cm)') );
    title('Total error (abs. val.) for (x,y,z) states')
    % Scale this plot to emphasize how small these errors are.
    %ylim([0 0.00125]);
    ylim([0 4]);
    % Make the font larger for these subplots that get squished.
    %set(gca,'FontSize',13);
    % Make the second plot:
    % Make a legend
    legend('No Dist.', 'With Dist.', 'Location', 'Northwest');
    hold off;
    subplot(2, 1, 2);
    hold on;
    % Adjust this error to degrees squared, since that's more intuitive than rad^2.
    % That means multiply by (180/pi)^2
    plot(t, sqrt(tracking_angle_squared_sum1)*180/pi, 'g.-');
    plot(t, sqrt(tracking_angle_squared_sum2)*180/pi, 'm.-');
    xlabel('Time (msec)');
    ylabel( sprintf('Sum Err. (deg)') );
    title('Total error (abs. val.) for angle (T, G, P) states');
    % Make the font larger for these subplots that get squished.
    %set(gca,'FontSize',13);
    % Scale:
    %ylim([0 0.05]);
    ylim([0 15]);
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






