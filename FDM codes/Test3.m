    
%---ASR (NODEs 2020)
% Shooting method Example 7.2 pg 180
% uses classicalRK_4 for IVPs
clc

global Re B0 Bb
tspan=[0 1]; % set time interval
Re=10;
Ba = [1, 0, 0, 0;
        0, 1, 0, 0;
        0, 0, 0, 0;
        0, 0, 0, 0];
    Bb = [0, 0, 0, 0;
        0, 0, 0, 0;
        1, 0, 0, 0;
        0, 1, 0, 0];


alpha = [0 0 1 0]';

I = eye(4);

figure1 = figure('Units','inches','Position',[0,0,6,6*0.7]); 

% first shooting
disp('first shooting')
disp('s is:')
s0  = [0 0 0 0]';
disp(s0)

[t,y1]=classicalRK_4(@combofun,tspan, [s0; I(:,1); I(:,2);I(:,3); I(:,4)], 100);
%[t,y2]=classicalRK_4(@combofun,tspan, [s0; I(:,2)], 100);
%[t,y3]=classicalRK_4(@combofun,tspan, [s0; I(:,3)], 100);


residual = B0*[y1(1:4,1)]+Bb*[y1(1:4,101)]-alpha;
disp('residual is:')
disp(residual)

p1 = plot(t,y1(1,:),'g-.');

zb = [y1(5:9,101), y1(10:14,101), y1(15:19,101), y1(20:24,101)];
s0 = s0-(B0+Bb*zb)\residual;
disp(s0)




function da = combofun(t,u)
Re=10;     
A = [0, 1, 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1;
            -Re*u(4), Re*u(3), Re*u(2), -Re*u(1)]; % Jacobian df/dy
dy = [u(2), u(3), u(4),Re*(u(2)*u(3)-u(1)*u(4))];
dz1 = A*[u(5:9)];
dz2 = A*[u(10:14)];
dz3 = A*[u(15:19)];
dz4= A*[u(20:24)];

da = [dy;dz1;dz2;dz3;dz4];
end
