close all; clear; clc;

% Set the print level to 5 in ipopt
% This script reads iterations, constraint violation, objective value
% from ipopt.out and visualizes them.

num_trails = 100;

data1 = zeros(num_trails,2001);
data2 = zeros(num_trails,2001);
data3 = zeros(num_trails,2001);

%%
for i = 1:num_trails
    ipopt_filename = sprintf('ipopt_digit-modified-Bezier-%d-N8.out', i);
    
    fid = fopen(ipopt_filename, 'r');
    
    iterations = [];
    constrViolation = [];
    objValue = [];
    
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
        end
    end
    
    fclose(fid);

    len = length(iterations);
    data1(i,1:len) = iterations;
    data2(i,1:len) = objValue;
    data3(i,1:len) = constrViolation;
end

%
figure; hold on;
for i = 1:num_trails
    plot(1:2001, log(data2(i,:)) / log(10), 'b');
end
xlim([1,2000]);
xlabel('iterations');
ylabel('log_{10} cost function');
title('cost function vs. iterations');

figure; hold on;
for i = 1:num_trails
    plot(1:2001, log(data3(i,:)) / log(10), 'b');
end
xlim([1,2000]);
xlabel('iterations');
ylabel('log_{10} constraint violation');
title('constraint violation vs. iterations');

figure;
plot(log(sort(data3(:,end))) / log(10))
xlabel('index of the experiment');
ylabel('log_{10} of the constraint violation');
title('final constraint violation for 100 experiments');
