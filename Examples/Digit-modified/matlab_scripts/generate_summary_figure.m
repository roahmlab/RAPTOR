clear; clc;

AA = load('comp-IDTO-Bezier-N8.mat');
BB = load('comp-IDTO-Bezier-N16.mat');
CC = load('comp-TROPIC.mat');

figure; hold on;
plot(sort(AA.data1));
plot(sort(BB.data1));
plot(sort(CC.data1));
legend('IDTO,Bezier,Chebshev,N=8', ...
       'IDTO,Bezier,Chebshev,N=16', ...
       'TROPIC');
xlabel('index of experiment');
ylabel('number of iterations to terminate');