% generate the reference traj used in ACC 2017

spacing = 0.1;
num_pts = 80;
direction = -1;

[traj_acc2017, ~] = get_ref_traj_invkin_constR_XZG(spacing, num_pts, direction);