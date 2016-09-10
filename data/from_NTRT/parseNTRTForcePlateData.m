function [ fpdata ] = parseNTRTForcePlateData( logfile_base, logfile_timestamp, calib_factor, make_plots )
%parseNTRTForcePlateData.m
%   Parses and plots the data from a single run of the
%   NASA Tensegrity Robotics Toolkit's AppHorizontalSpine with 4 x ForcePlateModel and Sensor
%   Andrew P. Sabelhaus
%   Berkeley Emergent Space Tensegrities Lab
%   Sept. 9, 2016
%
%   @param[in] logfile_base, a string. Same string as is passed in to ForcePlateSensor in NTRT.
%       For the ICRA 2017 paper, this will be ($NTRT_BASE)/resources/src/forcePlate/forcePlateDemo/logs/
%   @param[in] logfile_timestamp, also a string, of the timestamp that's used in the name of the logfiles.
%       See the log files themselves for information about the timestamp structure: it's both the date and time.
%   @param[in] calib_factor, the calibration factor to use for the force data. For example, if 
%       NTRT is set to gravity = 98.1, and the structure is scaled to scale = 0.5, then
%       calib_factor should be 10 * 0.5 = 5.
%   @param[in] make_plots, a flag that controls creation of graphs of the data or not
%   @retvar[out] fpdata, a cell array will all the force plate data nicely organized.

% First, construct the paths to the four force plates,
% hard-coding the labels I used for each of those plates.
% Cell arrays are nice here.
% Order is: RearLeft, RearRight, FrontLeft, FrontRight, TotalForces
fpdata = {};
fpdata{1}.path = strcat( logfile_base, 'ForcePlateSensor_FP_RearLeft_', ...
    logfile_timestamp, '.txt');
fpdata{2}.path = strcat( logfile_base, 'ForcePlateSensor_FP_RearRight_', ...
    logfile_timestamp, '.txt');
fpdata{3}.path = strcat( logfile_base, 'ForcePlateSensor_FP_FrontLeft_', ...
    logfile_timestamp, '.txt');
fpdata{4}.path = strcat( logfile_base, 'ForcePlateSensor_FP_FrontRight_', ...
    logfile_timestamp, '.txt');

% The dimensions to use in csvread: 
% Four columns, Time | Fx | Fy | Fz
% Data starts at the second row (row one.)

% For each of these, read in the data for each of these logs
for i=1:4
    % Call csvread
    fpdata{i}.data = csvread( fpdata{i}.path, 1, 0);
    % Adjust the data. All rows except the timestamp
    fpdata{i}.data(:, 2:end) = fpdata{i}.data(:, 2:end) .* (1/calib_factor);
end

% Combine the results, maybe this is used for testing of the whole robot
num_samples = size(fpdata{1}.data, 1);
% Record the totals for x, y, and z directions, plus time.
fpdata{5}.data = zeros(num_samples, 4);
% But we also need to copy over the timestamps for column 1.
% Arbitrarily choose the first set of data from which to copy the timestamps.
fpdata{5}.data(:,1) = fpdata{1}.data(:,1);
for i=1:4
    % Copy over the x, y, and z, columns 2, 3, 4
    fpdata{5}.data(:,2) = fpdata{5}.data(:,2) + fpdata{i}.data(:,2);
    fpdata{5}.data(:,3) = fpdata{5}.data(:,3) + fpdata{i}.data(:,3);
    fpdata{5}.data(:,4) = fpdata{5}.data(:,4) + fpdata{i}.data(:,4);
end

if( make_plots )
    % Plot the forces in x, y, z for each forceplate.
    
    % X:
    figure;
    hold on;
    for i=1:4
        % plot t vs. x for each plate
        plot( fpdata{i}.data(:,1), fpdata{i}.data(:,2) )
    end
    title('Force plate Fx forces vs. time');
    ylabel('Force Fx (N)');
    xlabel('Time (sec)');
    legend('RearLeft', 'RearRight', 'FrontLeft', 'FrontRight' );
    hold off;
    
    % Y
    figure;
    hold on;
    for i=1:4
        % plot t vs. y for each plate
        plot( fpdata{i}.data(:,1), fpdata{i}.data(:,3) )
    end
    title('Force plate Fy forces vs. time');
    ylabel('Force Fy (N)');
    xlabel('Time (sec)');
    legend('RearLeft', 'RearRight', 'FrontLeft', 'FrontRight' );
    hold off;
    
    % Z
    figure;
    hold on;
    for i=1:4
        % plot t vs. z for each plate
        plot( fpdata{i}.data(:,1), fpdata{i}.data(:,4) )
    end
    title('Force plate Fz forces vs. time');
    ylabel('Force Fz (N)');
    xlabel('Time (sec)');
    legend('RearLeft', 'RearRight', 'FrontLeft', 'FrontRight' );
    hold off;
    
    % Plot the total Fy, for perspective on how much the robot weighs.
    figure;
    hold on;
    plot( fpdata{5}.data(:,1), fpdata{5}.data(:,3) );
    title('Total forces in Y, NTRT force plates, vs. time');
    ylabel('Force Fy (N)');
    xlabel('Time (sec)');
    hold off;
end




