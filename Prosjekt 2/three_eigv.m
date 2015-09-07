clear all

% This program plots the probability density |u(r)|^2, where u is the
% wavefunction

fileID = fopen('eig1.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
eig1 = cols(2,:);
x = cols(1,:); % Just need to do this once
fclose(fileID);

fileID = fopen('eig2.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
eig2= cols(2,:);
fclose(fileID);

fileID = fopen('eig3.dat','r');
formatSpec = '%f %f';
size = [2 Inf];
cols = fscanf(fileID,formatSpec,size);
eig3 = cols(2,:);
fclose(fileID);


Figure1=figure(1);
set(Figure1,'defaulttextinterpreter','latex');
hold all
plot(x,eig1.^2)
plot(x,eig2.^2)
plot(x,eig3.^2)
title('Three lowest energy states, $|eigvec_n(\rho)|^2$, n=1,2,3','FontSize',14)
xlabel('$\rho$','FontSize',13)
ylabel('$|eigv_n(\rho)|^2$, n=1,2,3','FontSize',13)
legend('n=1, Lowest energy state, n=1','n=2, Second lowest energy state, n=2','n=3, Third lowest energy state')
