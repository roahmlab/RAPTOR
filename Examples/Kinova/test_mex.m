clear; clc;

q_0 = zeros(7,1);
q_dot_0 = zeros(7,1);
q_ddot_0 = zeros(7,1);
z0 = zeros(7,1); % initial guess of k center
duration = 1.4;
q_des = ones(7,1);
obstacles{1} = zonotope(zeros(3,1), 0.1 * eye(3));
obstacles{2} = zonotope([0;0;1.8], 0.1 * eye(3));
joint_limits_buffer = 0.01 * ones(7,1);
torque_limits_buffer = 10 * ones(7,1);

cpp_input_file = fopen('../../build/oracle_input_buffer.txt', 'w');

for ind = 1:length(q_0)
    fprintf(cpp_input_file, '%.10f ', q_0(ind));
end
fprintf(cpp_input_file, '\n');
for ind = 1:length(q_dot_0)
    fprintf(cpp_input_file, '%.10f ', q_dot_0(ind));
end
fprintf(cpp_input_file, '\n');
for ind = 1:length(q_ddot_0)
    fprintf(cpp_input_file, '%.10f ', q_ddot_0(ind));
end
fprintf(cpp_input_file, '\n');
fprintf(cpp_input_file, '%.10f\n', duration);
for ind = 1:length(z0)
    fprintf(cpp_input_file, '%.10f ', z0(ind)); % initial guess of k center
end
fprintf(cpp_input_file, '\n');
for ind = 1:length(q_des)
    fprintf(cpp_input_file, '%.10f ', q_des(ind));
end
fprintf(cpp_input_file, '\n');
fprintf(cpp_input_file, '%d\n', max(length(obstacles), 0));
for obs_ind = 1:length(obstacles)
    temp = reshape(obstacles{obs_ind}.Z, [1,size(obstacles{obs_ind}.Z,1) * size(obstacles{obs_ind}.Z,2)]);
    for ind = 1:length(temp)
        fprintf(cpp_input_file, '%.10f ', temp(ind));
    end
    fprintf(cpp_input_file, '\n');
end
for ind = 1:length(joint_limits_buffer)
    fprintf(cpp_input_file, '%.10f ', joint_limits_buffer(ind));
end
fprintf(cpp_input_file, '\n');
for ind = 1:length(torque_limits_buffer)
    fprintf(cpp_input_file, '%.10f ', torque_limits_buffer(ind));
end
fprintf(cpp_input_file, '\n');

fclose(cpp_input_file);

% call cuda program in terminal
% you have to be in the proper path!
terminal_output = system('./../../build/Kinova_example_mex');

if terminal_output == 0
    k_center_opt = readmatrix('solution-kinova.txt', 'FileType', 'text');

    if length(k_center_opt) == 1
        disp('Unable to find new trajectory!');
    else
        disp('New trajectory found!');
    end
else
    error('CUDA program error! Check the executable path in armour-dev/kinova_src/kinova_simulator_interfaces/uarmtd_planner');
end