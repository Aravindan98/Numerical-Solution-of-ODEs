ti=ones(2,1)';
vi=zeros(2,1)';

Ba=zeros(4);
Ba(1,1)=(2/(pi*ti(1)))^(1/2);
Ba(1,2)=1;
Ba(1,3)=-(vi(1)*(1/ti(1))^(3/2))/sqrt(2*pi);
Ba(3,3)=(2/pi)^(1/2)*(1/ti(1))^(3/2)*(ti(1)+1);
Ba(3,4)=1;

Bb=zeros(4);
Bb(2,1)=-(2/(pi*ti(2)))^(1/2);
Bb(2,2)=1;
Bb(4,3)=-(2/pi)^(1/2)*(1/ti(2))^(3/2)*(ti(2)+1);
Bb(2,3)=(vi(2)*(1/ti(2))^(3/2))/sqrt(2*pi);
Bb(4,4)=1;

% ASR NODES 2020
% The non-linear mid-point method
% SOLVES y'(x)=q(x,y) 
% g(y(a),g(b)) = 0 
%Ex 6
%EOC
format long;
close all;
clear all;

Nf = 800;
Kn=[0.0]
[xf, y1f, y2f, y3f, y4f] = ass2(Nf,0.068); 
mesh = [25,50,100,200]';
ratio =2*[16,8,4,2]';
Error = zeros(4,4);


for i = 1:4
    [x, y1, y2, y3, y4] = ass2(mesh(i),0.068);
    k = ratio(i);
    disp([size(y1),size(y1f)])
    Error(i,1) = max(abs(y1f(1:k:Nf+1)-y1));
    Error(i,2) = max(abs(y2f(1:k:Nf+1)-y2));
    Error(i,3) = max(abs(y3f(1:k:Nf+1)-y3));
    Error(i,4) = max(abs(y4f(1:k:Nf+1)-y4));
end

p1 = polyfit(log(1./mesh),log(Error(:,1)),1);
p2 = polyfit(log(1./mesh),log(Error(:,2)),1);
p3 = polyfit(log(1./mesh),log(Error(:,3)),1);
p4 = polyfit(log(1./mesh),log(Error(:,4)),1);


figure(1);
loglog(mesh,Error(:,1),'-ok','LineWidth',1)
hold on;
loglog(mesh,Error(:,2),'-sr','LineWidth',1);
loglog(mesh,Error(:,3),'-^b','LineWidth',1)
loglog(mesh,Error(:,4),'-*m','LineWidth',1)
grid on;

xlabel('$N$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$\epsilon_{h}$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
