clear all

% This program plots the probability density |u(r)|^2, where u is the
% wavefunction

fileID = fopen('rep_w0_01.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
rep_w0_01 = cols(2,:);
x = cols(1,:); % Just need to do this once
fclose(fileID);

fileID = fopen('rep_w0_5.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
rep_w0_5 = cols(2,:);
fclose(fileID);

fileID = fopen('rep_w1.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
rep_w1 = cols(2,:);
fclose(fileID);

fileID = fopen('rep_w5.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
rep_w5 = cols(2,:);
fclose(fileID);

fileID = fopen('norep_w0_01.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
w0_01 = cols(2,:);
fclose(fileID);

fileID = fopen('norep_w0_5.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
w0_5 = cols(2,:);
fclose(fileID);

fileID = fopen('norep_w1.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
w1 = cols(2,:);
fclose(fileID);

fileID = fopen('norep_w5.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
w5 = cols(2,:);
fclose(fileID);

Figure1=figure(1);
set(Figure1,'defaulttextinterpreter','latex');

subplot(4,1,1)
plot(x,rep_w0_01.^2,'b')
hold on
plot(x,w0_01.^2,'r')
title('Probability density for different $\omega_r$ with and w/o repulsion','FontSize',14)
xlabel('$\rho$','FontSize',13)
ylabel('$|u(\rho)|^2$','FontSize',13)
legend('With repulsion, \omega_r=0.01','No repulsion, \omega_r=0.01')

subplot(4,1,2)
plot(x,rep_w0_5.^2,'b')
hold on
plot(x,w0_5.^2,'r')
xlabel('$\rho$','FontSize',13)
ylabel('$|u(\rho)|^2$','FontSize',13)
legend('With repulsion, \omega_r=0.5','No repulsion, \omega_r=0.5')

subplot(4,1,3)
plot(x,rep_w1.^2,'b')
hold on
plot(x,w1.^2,'r')
xlabel('$\rho$','FontSize',13)
ylabel('$|u(\rho)|^2$','FontSize',13)
legend('With repulsion, \omega_r=1','No repulsion, \omega_r=1')

subplot(4,1,4)
plot(x,rep_w5.^2,'b')
hold on
plot(x,w5.^2,'r')
xlabel('$\rho$','FontSize',13)
ylabel('$|u(\rho)|^2$','FontSize',13)
legend('With repulsion, \omega_r=5','No repulsion, \omega_r=5')






