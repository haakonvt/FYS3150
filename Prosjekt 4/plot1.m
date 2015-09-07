clear all
clc

% Plot of spatial concentration neurotransmitters after different time

Nx = 10;
Nt = 200;

x  = linspace(0,1,Nx);
xx = linspace(0,1,Nx*10);

fileID = fopen('ExplForwEuler.dat','r');
sol_EFE = fscanf(fileID,'%f %f',[Nt Nx])';
fclose(fileID);

fileID = fopen('ImplBackEuler.dat','r');
sol_IBE = fscanf(fileID,'%f %f',[Nt Nx])';
fclose(fileID);

fileID = fopen('u_analytic.dat','r');
sol_exact = fscanf(fileID,'%f %f',[Nt Nx*10])';
fclose(fileID);

subplot(3,1,1)
Figure1=figure(1);
set(Figure1,'defaulttextinterpreter','latex');
hold all
for i =[linspace(2,20,10) 30 50 100 150 200]
    plot(x,sol_EFE(:,i))
end
%legend('Forw, t = 0.02','Forw, t = 0.05','Forw, t = 0.1','Forw, t = 0.5')
title('Explicit Forward Euler','FontSize',15)
xlabel('x','FontSize',13)
ylabel('Concentration','FontSize',13)
axis([0 1 0 1])



subplot(3,1,2)
hold all
for i =[linspace(2,20,10) 30 50 100 150 200]
    plot(x,sol_IBE(:,i))
end

title('Implicit Back Euler','FontSize',15)
%legend('Back, t = 0.02','Back,t = 0.05','Back,t = 0.1','Back,t = 0.5')
xlabel('x','FontSize',13)
ylabel('Concentration','FontSize',13)
axis([0 1 0 1])



subplot(3,1,3)
hold all
for i =[linspace(2,20,10) 30 50 100 150 200]
    plot(xx,sol_exact(:,i))
end

title('Analytic Solution','FontSize',15)
%legend('Back, t = 0.02','Back,t = 0.05','Back,t = 0.1','Back,t = 0.5')
xlabel('x','FontSize',13)
ylabel('Concentration','FontSize',13)
axis([0 1 0 1])








