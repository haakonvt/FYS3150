clear all

% This program need main.cpp to be run three times with diferent values of
% n, preferably n=10,100,1000

fileID = fopen('rel_error_n10.dat','r');
formatSpec = '%f';
N10 = fscanf(fileID,formatSpec);
fclose(fileID);

fileID = fopen('rel_error_n100.dat','r');
formatSpec = '%f';
N100 = fscanf(fileID,formatSpec);
fclose(fileID);

fileID = fopen('rel_error_n1000.dat','r');
formatSpec = '%f';
N1000 = fscanf(fileID,formatSpec);
fclose(fileID);

% h = 1.0/(n+1);
x10   = linspace(1/11,1-1/11,10);
x100  = linspace(1/101,1-1/101,100);
x1000 = linspace(1/1001,1-1/1001,1000);

figure(1)
subplot(3,1,3)
plot(x1000,N1000)
legend('Relative error for n=1000')
xlabel('x')
ylabel('log10((v-u)/u)')

% ---------------------------------

subplot(3,1,2)
plot(x100,N100)
legend('Relative error for n=100')
xlabel('x')
ylabel('log10((v-u)/u)')

% ---------------------------------

subplot(3,1,1)
plot(x10,N10)
legend('Relative error for n=10')
xlabel('x')
ylabel('log10((v-u)/u)')
title('Relative error for increasing n','FontSize',12)