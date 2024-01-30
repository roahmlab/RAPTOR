close all; clear; clc;

% Set the print level to 5 in ipopt
% This script reads iterations, constraint violation, objective value
% from ipopt.out and visualizes them.

num_trails = 100;

data1 = zeros(num_trails,1);
data2 = zeros(num_trails,2001);
data3 = zeros(num_trails,2001);
data4 = zeros(num_trails,1);
data5 = zeros(num_trails,1);
running_time = zeros(num_trails, 1);

%%
for i = 1:num_trails
%     ipopt_filename = sprintf('ipopt_digit-modified-Bezier-%d-N16-20-Chebyshev.out', i);
    ipopt_filename = sprintf('/home/roahmlab/Documents/TROPIC/examples/digit-modified/data/ipopt-%d-d20.out', i);

    fid = fopen(ipopt_filename, 'r');
    
    iterations = [];
    constrViolation = [];
    objValue = [];
    energy = 0;
    
    while ~feof(fid)
        line = fgetl(fid);
        line_process = strsplit(line, ' ');

        if isempty(line_process{1})
            line_process = line_process(2:end);
        end

        if length(line_process) == 10
            if line_process{1}(end) == 'r'
                iterations_text = line_process{1}(1:end-1);
            else
                iterations_text = line_process{1};
            end

            iteration_number = str2double(iterations_text);

            if isnan(iteration_number)
                continue;
            end

            iterations = [iterations, iteration_number];
            objValue = [objValue, str2double(line_process{2})];
            constrViolation = [constrViolation, str2double(line_process{3})];
        elseif length(line_process) == 6 && ~isempty(str2double(line_process{end})) && strcmp(line_process{2}, 'seconds')
            running_time(i) = str2double(line_process{end});
        elseif length(line_process) == 4 && strcmp(line_process{1}, 'Objective')
            energy = str2double(line_process{end});
        end
    end
    
    fclose(fid);

    len = length(iterations);
%     data1(i,1:len) = iterations;
    data1(i) = len;
    data2(i,1:len) = objValue;
    data3(i,1:len) = constrViolation;
    data4(i) = constrViolation(end);
    data5(i) = energy;
end

%%
% figure; hold on;
% for i = 1:num_trails
%     plot(1:2001, log(data2(i,:)) / log(10), 'b');
% end
% xlim([1,2000]);
% xlabel('iterations');
% ylabel('log_{10} cost function');
% title('cost function vs. iterations (TROPIC)');
% 
% figure; hold on;
% for i = 1:num_trails
%     plot(1:2001, log(data3(i,:)) / log(10), 'b');
% end
% xlim([1,2000]);
% xlabel('iterations');
% ylabel('log_{10} constraint violation');
% title('constraint violation vs. iterations (TROPIC)');
% 
% figure;
% plot(log(sort(data3(:,end))) / log(10))
% xlabel('index of the experiment');
% ylabel('log_{10} of the constraint violation');
% title('final constraint violation for 100 experiments (TROPIC)');

%%
% save('comp-Bezier-N16-50-Uniform.mat', 'data1', "data2", "data3", "data4");

disp(vpa([median(data4), max(data4)], 6));
disp([median(data1), max(data1)])
disp(vpa([median(data5), max(data5)], 6));
disp(vpa([median(running_time), max(running_time)], 6));
