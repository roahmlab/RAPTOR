close all; clear; clc;

robot = importrobot('digit-v3-modified.urdf');
robot.DataFormat = 'col';

traj = readmatrix('trajectory-digit-modified-Bezier.txt');

q = traj(1:20, :);
qd = traj(21:40, :);
qdd = traj(41:60, :);
u = traj(61:72, :);
lambda = traj(73:end, :);

for i = 1:size(q,2)
    show(robot, q(:,i), 'frames', 'off', 'Collisions','off');
    hold on;
end
axis([-1 1 -1 1 0 2]);