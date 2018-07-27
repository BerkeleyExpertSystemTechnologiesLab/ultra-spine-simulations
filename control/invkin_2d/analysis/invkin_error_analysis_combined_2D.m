% invkin_error_analysis_combined_2D.m
% Copyright 2018 Andrew P. Sabelhaus, Berkeley Emergent Space Tensegrities Lab
% This function loads in a set of saved data from multiple MPC runs and
% plots a comparison.
% For the 2D MPC results.

% July 2018 - for inverse kinematics.

function [errors1, errors2] = invkin_error_analysis_combined_2D( file_name1, file_name2, path_to_data_folder, plots_flag )
% Inputs:
%   file_name1,2 = name of the data file, needs to include '.mat'
%   path_to_data_folder = location of the data file to read in.
%   plots_flag = create plots (1) or do not make plots (0).
% Outputs:
%   none, since this function just plots.

% a bit of cleanup.
close all;

% Some hard-coded variables 
% TO-DO: pass these in to the function?
% Size of the text in the figures
fontsize = 12;

% Call mpc_error_analysis)2D to calculate all the errors for both MPC runs.
errors1 = invkin_error_analysis_2D( file_name1, path_to_data_folder);
errors2 = invkin_error_analysis_2D( file_name2, path_to_data_folder);

% Make plots, if flag is set
if plots_flag
    % Calculate the time in seconds at each point.
    % Assume that the times for errors1 are the same as for errors2.
    t = errors1.dt: errors1.dt : (errors1.num_pts-1)*errors1.dt;
    % Adjust time so that everything is in milliseconds. This makes things clearer.
    % Note that for the 2D model, we're looking at dt=1e-5, so even though
    % we're using millisec, will be fractions of msec vs. just msec in
    % comparison to 3D.
    t = t * 1000;
    % A good figure window setup is 'Position',[100,100,500,300].
    size(t)
    
    % Pick out some of the variables so that it's easier to write the code below:
    % Assume that both the data files use the same reference trajectory.
    ref_traj = errors1.ref_traj;
    result_traj1 = errors1.result_traj;
    result_traj2 = errors2.result_traj;
    
    %% Plot the XZ positions of the reference trajectory and resulting trajectory.
    % This should give a good visual representation of the spine's movement.
    
    % Single vertebra:
    trajectories_handle = figure;
    hold on;
    set(trajectories_handle,'Position',[100,100,500,300]);
    % Plot the reference: X vs. Z. These are states 1 and 2.
    % Scale the lengths here to get cm.
    plot( ref_traj(1,:)*100, ref_traj(2, :)*100, 'b','LineWidth',2);
    % Plot the output of the MPC:
    % Scale the lengths here to get cm.
    % 2018-07-27: switched ordering to put noise-less version "in front"
    plot( result_traj2(1,:)*100, result_traj2(2,:)*100, 'm','LineWidth',2);
    plot( result_traj1(1,:)*100, result_traj1(2,:)*100, 'g','LineWidth',2);
    legend('Target Trajectory', 'Result, With Noise', 'Result, No Noise', 'Location','Southeast');
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of Vertebra, Inv. Kin. Control');
    set(gca,'FontSize',fontsize);
    % Set the limits more intelligently.
    xlim([-6.05, 0]);
    %ylim([29.1, 30.05]);
    % for a better sense of how big the model error is:
    ylim([27.5, 30.5]);
    % Scale the plot?
    hold off;
    
    %% Plot the individual errors for X, Z, gamma
    
    % Pick out the tracking errors.
    tracking_error1 = errors1.tracking_error;
    tracking_error2 = errors2.tracking_error;
    
    % Create the handle for the overall figure
    % Also, use the OpenGL renderer so that symbols are formatted correctly.
    %errors_handle = figure('Renderer', 'opengl');
    errors_handle = figure;
    hold on;
    set(gca,'FontSize',fontsize);
    % This figure will have 3 smaller plots, so make it larger than my
    % usual window dimensions.
    % Was:
    %set(errors_handle,'Position',[100,100,900,300]);
    % So, just halve the vertical coordinate., so... 450?
    set(errors_handle,'Position',[100,100,450,300]);
    % Start the first subplot
    subplot_handle = subplot(3, 1, 1);
    hold on;
    
    % Plot the X errors
    % Scale the lengths here to get centimeters.
    % Ignore the first datapoint (meaningless)
    plot(t, tracking_error1( 1, 2:end)*100, 'Color', 'g', 'LineWidth', 2);
    % for the disturbances, too
    plot(t, tracking_error2( 1, 2:end)*100, 'Color', 'm', 'LineWidth', 2);
    % Plot the zero line
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %xlabel('Time (msec)');
    ylabel('e_X (cm)');
    % Only create a title for the first plot, that will serve for all the others too.
    %title('Tracking Errors in X Y Z  \theta \gamma \psi');
    title('   State Errors, Input Reference Tracking');
    % Scale the plot. A good scale here is...
    %ylim([-2.0 3.5]);
    % Adjust by roughly the amount we scaled the disturbances: 1/6 of the
    % length. Plus a small change to make the numbers prettier, about -(1/3)+0.05
    ylim([-0.1, 0.55]);
    
    % Make the legend
    nodisturblabel = sprintf('No disturb.');
    disturblabel = sprintf('With disturb.');
    legend_handle = legend(nodisturblabel, disturblabel, 'Location', 'North', 'Orientation', 'horizontal');
    
    % Move the plot very slightly to the left
    %P = get(subplot_handle,'Position')
    %set(subplot_handle,'Position',[P(1)-0.05 P(2) P(3) P(4)])
    %set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.05 P(3)+0.01 P(4)-0.04])
    hold off;
    
    % Title the whole plot
    %mtit('Tracking Errors, With Disturbances','xoff',0.5,'yoff',0.2)
    
    % move the legend to the corner of the figure.
    % Position is: [left, bottom, width, height]
    % These units are expressed as percents of the total figure
    % Bottom left, for the 1x6 plot:
    %legend_position = [0.22 0.042 0 0];
    % Bottom right, 1x6 plot:
    %legend_position = [0.82 0.042 0 0];
    % Bottom right, 3x2 plot:
    %legend_position = [0.935 0.6 0.001 0];
    %set(legend_handle,'Position', legend_position);
    %title('Tracking Error in \psi');
    

    % Plot the Z errors
    subplot_handle = subplot(3,1,2);
    hold on;
    % Ignore the first datapoint (meaningless)
    plot(t, tracking_error1( 2, 2:end)*100, 'Color', 'g', 'LineWidth', 2);
    % for the disturbances, too
    plot(t, tracking_error2( 2, 2:end)*100, 'Color', 'm', 'LineWidth', 2);
    % Plot the zero line
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    %xlabel('Time (msec)');
    ylabel('e_Z (cm)');
    %title('Tracking Error in Y');
    % Adjust by roughly the amount we scaled the disturbances: 1/6 of the
    % length. Plus a small change to make the numbers prettier, about -(1/3)+0.05
    %ylim([-0.2, (7/12)-0.1]); 
    ylim([-0.1, 0.55]);

    % Move the plot very slightly to the left
    % For these lower figures, move them upwards a bit more.
    %P = get(subplot_handle,'Position')
    %set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.07 P(3)+0.01 P(4)-0.04])
    
    % Move it slightly down.
    %legend_position = [0.935 0.6 0.001 0];
    %set(legend_handle,'Position', legend_position);
    
    
    hold off;
    
    % Plot the gamma errors
    subplot_handle = subplot(3,1,3);
    hold on;
    % Ignore the first datapoint (meaningless)
    % Convert rad to degrees
    plot(t, tracking_error1( 3, 2:end) * (180/pi), 'Color', 'g', 'LineWidth', 2);
    % for the disturbances, too
    plot(t, tracking_error2( 3, 2:end) * (180/pi), 'Color', 'm', 'LineWidth', 2);
    % Plot the zero line
    %plot(t, zero_line, 'b-', 'LineWidth','1');
    refline_handle = refline(0,0);
    set(refline_handle, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
    %legend('Vertebra 1 (Bottom)', 'Vertebra 2 (Middle)', 'Vertebra 3 (Top)');
    %xlabel('Time (msec)');
    ylabel('e_\gamma (deg)');
    %title('Tracking Error in Y');
    % Angles limits: maybe the same as in 3D?
    ylim([-5 1]);
    % Move the plot very slightly to the left
    % For these lower figures, move them upwards a bit more.
    %P = get(subplot_handle,'Position')
    %set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.07 P(3)+0.01 P(4)-0.04])
    
    % Finally, a label in X at the bottom
    xlabel('Time (msec)');
    

    hold off;
    
     
    % For Z: Move the plot very slightly to the left
    %P = get(subplot_handle,'Position')
    %set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.10 P(3)+0.01 P(4)-0.04])
    

    % for the right-hand side ones, was
    %P = get(subplot_handle,'Position')
    %set(subplot_handle,'Position',[P(1)-0.06 P(2)+0.05 P(3)+0.01 P(4)-0.04])
    


end

end






