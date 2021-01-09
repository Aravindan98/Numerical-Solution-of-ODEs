% ASR NODES 2020
% The non-linear mid-point method
% SOLVES y'(x)=q(x,y) 
% g(y(a),g(b)) = 0 
%Ex 6
function [x, y1, y2, y3, y4] = ass2(N)
clear all
format long
clc

N = 200; % no of discritization points
mm = 4;
k = 0.1;
F = 0.23;
h = 1/N; % mesh size
%Re = 10;%Reynolds number

x = (0:h:1)'; % nodes including boundary nodes
M = zeros(mm*(N+1),mm*(N+1));
b = zeros(mm*(N+1),1);

Done = 0; 
Max_iter = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0 = zeros(mm*(N+1),1);
for i = 1:N+1
   y0(mm*(i-1)+1) = 0;
   y0(mm*(i-1)+2) = 0;
   y0(mm*(i-1)+3) = 1;
   y0(mm*(i-1)+4) = 0;
end

iter = 0;

Ba = [(2/(pi*y0(3)))^(1/2), 1, -(y0(1)*(1/y0(3))^(3/2))/sqrt(2*pi), 0;
        0, 0, 0, 0;
        0, 0, (2/pi)^(1/2)*(1/y0(3))^(3/2)*(y0(3)+1), 0;
        0, 0, 0, 0];
    Bb = [0, 0, 0, 0;
        -(2/(pi*y0(mm*N+3)))^(1/2), 1, (y0(mm*N+1))*(1/y0((mm*N+1))^(3/2))/sqrt(2*pi), 0;
        0, 0, 0, 0;
        0, 0, -(2/pi)^(1/2)*((1/y0(mm*N+3))^(3/2))*(y0(mm*N+3)+1), 1];

while (~Done)
    for i = 1:N
		
        xip1 = 0.5*(x(i)+x(i+1));
        yip1 = 0.5*(y0(mm*(i-1)+1:mm*i)+y0(mm*i+1:mm*(i+1)));
		Aip1 = [0, -1/(k), 0, 0;
            0, 0, -F/(yip1(3)*yip1(3)), 0;
            0, 0, 0, -4/(15*k);
            0, 2*yip1(2)/(k), 0, 0];
        
        Qip1 = [-yip1(2)/(k); F/yip1(3); -4*yip1(4)/(15*k); yip1(2)*yip1(2)/(k)];
        fip1 = Qip1-Aip1*yip1;
        M(mm*(i-1)+1:mm*i,mm*(i-1)+1:mm*i) = -(0.5*Aip1+eye(mm)/h);
        M(mm*(i-1)+1:mm*i,mm*i+1:mm*(i+1)) = (-0.5*Aip1+eye(mm)/h);
        b(mm*(i-1)+1:mm*i) = fip1;

    end 
		[Ba,Bb,alpha]=getboundary(@g,y0,N,Ba,Bb,mm);
         M(mm*((N+1)-1)+1:mm*(N+1),mm*((N+1)-1)+1:mm*(N+1)) = Bb;
		M(mm*((N+1)-1)+1:mm*(N+1),1:mm) = Ba;
		b(mm*((N+1)-1)+1:mm*(N+1)) = alpha;
       
    y = M\b;
    
    error = max(abs(y-y0));% max |y_i-y_i^*|
    iter = iter+1;
    disp(['Error in step ',num2str(iter),' is : ',num2str(error)])
    y0 = y;

    Done = (iter>=Max_iter)||(error<0.000001);        
end


y1 = zeros(N+1,1);
y2 = zeros(N+1,1);
y3 = zeros(N+1,1);
y4 = zeros(N+1,1);
for i = 1:N+1
    y1(i) = y(mm*(i-1)+1);
    y2(i) = y(mm*(i-1)+2);
    y3(i) = y(mm*(i-1)+3);
    y4(i) = y(mm*(i-1)+4);
end

figure(1);
plot(x,y1,'-b','LineWidth',1),grid on;
xlabel('$x$','FontSize',13,'Color','k', 'Interpreter', 'latex');
ylabel('$u(x)$','FontSize',13,'FontWeight','bold','Color','b', 'Interpreter', 'latex');

figure(2);
plot(x,y2,'-r','LineWidth',1),grid on;
xlabel('$x$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$u(x)$','FontSize',13,'FontWeight','bold','Color','b', 'Interpreter', 'latex')

figure(3);
plot(x,y3,'-g','LineWidth',1),grid on;
xlabel('$x$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$u(x)$','FontSize',13,'FontWeight','bold','Color','b', 'Interpreter', 'latex')

figure(4);
plot(x,y4,'-k','LineWidth',1),grid on;
xlabel('$x$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$u(x)$','FontSize',13,'FontWeight','bold','Color','b', 'Interpreter', 'latex')
end

function [Ba,Bb,alpha]=getboundary(g,y,N,Ba,Bb,mm)
	ti=[y(3),y(mm*N+3)]';
vi=[y(1),y(mm*N+1)]';
Ba(1,2)=1;
Ba(1,3)=-(vi(1)*(1/ti(1))^(3/2))/sqrt(2*pi);
Ba(3,3)=(2/pi)^(1/2)*(1/ti(1))^(3/2)*(ti(1)+1);
Ba(3,4)=1;

Bb(2,1)=-(2/(pi*ti(2)))^(1/2);
Bb(2,2)=1;
Bb(4,3)=-(2/pi)^(1/2)*(1/ti(2))^(3/2)*(ti(2)+1);
Bb(2,3)=(vi(2)*(1/ti(2))^(3/2))/sqrt(2*pi);
Bb(4,4)=1;

alpha = Ba*[y(1),y(2),y(3),y(4)]' + Bb*[y(mm*N+1),y(mm*N+2),y(mm*N+3),y(mm*N+4)]' - feval(g,y,mm,N);
end

function [res] = g(y,mm,N)
res = zeros(mm,1);
res(1) = y(2) + ((2/(pi*y(3)))^(1/2))*y(1);
res(2) = y(mm*N+2) - ((2/(pi*y(mm*N+3)))^(1/2))*y(mm*N+1);
res(3) = y(4) + 2*((2/(pi*y(3)))^(1/2))*(y(3)-1);
res(4) = y(mm*N+4) - 2*((2/(pi*y(mm*N+3)))^(1/2))*(y(mm*N+3)-1);
end
