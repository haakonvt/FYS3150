clear all
pause on

% This program plots the orbits of one or more planets

year = 100;
N    = 120*year;

fileID = fopen('x.dat','r');
x = fscanf(fileID,'%f %f',[N 10])';
fclose(fileID);

fileID = fopen('y.dat','r');
y = fscanf(fileID,'%f %f',[N 10])';
fclose(fileID);

fileID = fopen('xx.dat','r');
xx = fscanf(fileID,'%f %f',[N 10])';
fclose(fileID);

fileID = fopen('yy.dat','r');
yy = fscanf(fileID,'%f %f',[N 10])';
fclose(fileID);

fileID = fopen('energy.dat','r');
E = fscanf(fileID,'%f %f',[N 3])';
fclose(fileID);

Figure1=figure(1);
set(Figure1,'defaulttextinterpreter','latex');

subplot(2,1,1)
hold all
for i=1:10
    plot(x(i,:),y(i,:))
end
legend('Sun, RK4', 'Earth, RK4','Jupiter (1000X mass) RK4')
%title('Sun-Earth-Jupiter (4 years, 5000 steps/year)','FontSize',14)
legend('Sun, RK4', 'Mercury, RK4','Venus, RK4','Earth, RK4','Mars, RK4','Jupiter, RK4','Saturn, RK4','Uranus, RK4','Neptune, RK4','Pluto, RK4')
axis equal
xlabel('X-axis, distance in $AU$','FontSize',13)
ylabel('Y-axis, distance in $AU$','FontSize',13)
title('Solar system','FontSize',14)

subplot(2,1,2)
hold all
for i=1:10
    plot(xx(i,:),yy(i,:))
end
%legend('Sun, VERLET', 'Earth, VERLET','Jupiter (1000X mass) VERLET')
legend('Sun, VER', 'Mercury, VER','Venus, VER','Earth, VER','Mars, VER','Jupiter, VER','Saturn, VER','Uranus, VER','Neptune, VER','Pluto, VER')
axis equal
xlabel('X-axis, distance in $AU$','FontSize',13)
ylabel('Y-axis, distance in $AU$','FontSize',13)

% ENERGY CONSERVATION
% Figure2=figure(2);
% set(Figure2,'defaulttextinterpreter','latex');
% NN = length(E(1,:))-5;
% time = linspace(0,year,NN);
% plot(time,E(1,1:NN),'b')
% hold on
% plot(time,E(2,1:NN),'r')
% hold on
% plot(time,E(1,1:NN)+E(2,1:NN),'k')
% hold on
% plot(time,E(3,1:NN),'g')
% %axis([0 100 -2.6E-3 -2.2E-3])
% legend('Potential energy','Kinetic energy','Total energy','Angular momentum')
% title('Conservation of energy and ang. momentum','FontSize',14)
% xlabel('Time [AU]','FontSize',13)
% ylabel('Pot. and kin.: [Joule], Ang. mom. [Joule$\cdot$ sec]','FontSize',13)


% ANIMATE SOLAR SYSTEM
% hold all
% axis equal
% counter = 3;
% when_draw = 3;
% pause(1)
% for i=N
%     if counter == when_draw;
%         plot(x(1,1:i),y(1,1:i),'kx')
%         plot(x(2,1:i),y(2,1:i),'b')
%         %plot(x(3,1:i),y(3,1:i),'r')
%         %plot(x(4,1:i),y(4,1:i),'k')
%         %plot(x(5,1:i),y(5,1:i),'m')
%         %plot(x(6,1:i),y(6,1:i),'b')
%         %plot(x(7,1:i),y(7,1:i),'r')
%         %plot(x(8,1:i),y(8,1:i),'k')
%         %plot(x(9,1:i),y(9,1:i),'m')
%         %plot(x(10,1:i),y(10,1:i),'b')
%         %axis([-2 40 -2 5])
%         drawnow
%         counter  = 0;
%     end
%     counter = counter + 1;
% end