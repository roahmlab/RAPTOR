close all; clear; clc;

AA = load('comp-Bezier-N8-50-Chebyshev.mat');
BB = load('comp-Bezier-N16-50-Chebyshev.mat');
CC = load('comp-Bezier-N8-50-Uniform.mat');
DD = load('comp-Bezier-N16-50-Uniform.mat');
EE = load('comp-TROPIC.mat');
FF = load('comp-crocoddyl-fwddyn.mat');
GG = load('comp-crocoddyl-invdyn.mat');

% figure; hold on;
% plot(log(sort(AA.data1));
% plot(log(sort(BB.data1));
% plot(log(sort(CC.data1));
% legend('IDTO,Bezier,Chebshev,N=8', ...
%        'IDTO,Bezier,Chebshev,N=16', ...
%        'TROPIC');
% xlabel('index of experiment');
% ylabel('number of iterations to terminate');

figure; hold on;
plot(log(sort(AA.data4(2:end))) / log(10));
plot(log(sort(BB.data4(2:end))) / log(10));
plot(log(sort(CC.data4(2:end))) / log(10));
plot(log(sort(DD.data4(2:end))) / log(10));
plot(log(sort(EE.data4(2:end))) / log(10));
plot(log(sort(FF.data3)) / log(10));
plot(log(sort(GG.data3)) / log(10));
ylim([-12 2]);
legend('IDTO,Bezier,Chebshev,N=8', ...
       'IDTO,Bezier,Chebshev,N=16', ...
       'IDTO,Bezier,Uniform,N=8', ...
       'IDTO,Bezier,Uniform,N=16', ...
       'TROPIC', ...
       'crocoddyl forward dynamics', ...
       'crocoddyl inverse dynamics');
xlabel('index of experiment');
ylabel('log_{10} constraint violation ratio');