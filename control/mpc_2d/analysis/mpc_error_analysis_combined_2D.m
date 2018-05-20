% mpc_error_analysis_combined_2D.m
% Copyright 2018 Andrew P. Sabelhaus, Berkeley Emergent Space Tensegrities Lab
% This function loads in a set of saved data from multiple MPC runs and
% plots a comparison.
% For the 2D MPC results.

function mpc_error_analysis_combined_2D( file_name1, file_name2, path_to_data_folder, plots_flag )
% Inputs:
%   file_name1,2 = name of the data file, needs to include '.mat'
%   path_to_data_folder = location of the data file to read in.
%   plots_flag = create plots (1) or do not make plots (0).
% Outputs:
%   none, since this function just plots.

% Some hard-coded variables 
% TO-DO: pass these in to the function?
% Size of the text in the figures
fontsize = 12;

% Call mpc_error_analysis)2D to calculate all the errors for both MPC runs.
errors1 = mpc_error_analysis_2D( file_name1, path_to_data_folder);
errors2 = mpc_error_analysis_2D( file_name2, path_to_data_folder);

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
    plot( result_traj1(1,:)*100, result_traj1(2,:)*100, 'g','LineWidth',2);
    plot( result_traj2(1,:)*100, result_traj2(2,:)*100, 'm','LineWidth',2);
    legend('Reference', 'Result, No Dist.', 'Result, With Dist.','Location','Southeast');
    xlabel('Position in X (cm)');
    ylabel('Position in Z (cm)');
    title('Position of Vertebra, Input Ref. Tracking');
    set(gca,'FontSize',fontsize);
    % Set the limits more intelligently.
    xlim([-2.05, 0]);
    ylim([9.79, 10.03]);
    % Scale the plot?
    hold off;

end

end






