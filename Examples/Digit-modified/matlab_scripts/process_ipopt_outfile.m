close all; clear; clc;

% Set the print level to 5 in ipopt
% This script reads iterations, constraint violation, objective value
% from ipopt.out and visualizes them.

data1 = zeros(100,1000);
data2 = zeros(100,1000);
data3 = zeros(100,1000);

for i = 1:100
    ipopt_filename = sprintf('ipopt_digit-modified-Bezier-%d.out', i);
    
    fid = fopen(ipopt_filename, 'r');
    
    iterations = [];
    constrViolation = [];
    objValue = [];
    
    while ~feof(fid)
        line = fgetl(fid);
        line_process = strsplit(line, ' ');
        if length(line_process) == 11
            if line_process{2}(end) == 'r'
                iterations = [iterations, str2double(line_process{2}(1:end-1))];
            else
                iterations = [iterations, str2double(line_process{2})];
            end
            objValue = [objValue, str2double(line_process{3})];
            constrViolation = [constrViolation, str2double(line_process{4})];
        end
    end
    
    fclose(fid);

    len = length(iterations);
    data1(i,1:len) = iterations;
    data2(i,1:len) = objValue;
    data3(i,1:len) = constrViolation;
end