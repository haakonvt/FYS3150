clear all

% This program need main.cpp to be run three times with diferent values of
% n, preferably n=10,100,1000

fileID = fopen('sol_vec_n10.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
N10 = fscanf(fileID,formatSpec,size);
fclose(fileID);

fileID = fopen('sol_vec_n100.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
N100 = fscanf(fileID,formatSpec,size);
fclose(fileID);

fileID = fopen('sol_vec_n1000.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
N1000 = fscanf(fileID,formatSpec,size);
fclose(fileID);

% Since values from fscanf is transposed, the rows must be used
x1000 = N1000(1,:); 
y = N1000(2,:);
yy = 1-(1-exp(-10)).*x1000 - exp(-10.*x1000);

figure(1)
subplot(3,1,3)
plot(x1000,y)
hold on
plot(x1000,yy,'r')
legend('Numeric solution for n=1000','Exact solution')
xlabel('x')
ylabel('-u(x)')
hold off

% ---------------------------------

x = N100(1,:); 
y = N100(2,:);

subplot(3,1,2)
plot(x,y)
hold on
plot(x1000,yy,'r')
legend('Numeric solution for n=100','Exact solution')
xlabel('x')
ylabel('-u(x)')
hold off

% ---------------------------------

x = N10(1,:); 
y = N10(2,:);

subplot(3,1,1)
plot(x,y)
hold on
plot(x1000,yy,'r')
legend('Numeric solution for n=10','Exact solution')
xlabel('x')
ylabel('-u(x)')
hold off
title('Convergence of numerical solution for increasing n')

