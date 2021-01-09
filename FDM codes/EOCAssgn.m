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
Kn=[0.068, 0.5];
[x1f, y11f, y12f, y13f, y14f] = ass2(Nf,Kn(1)); 
[x2f, y21f, y22f, y23f, y24f] = ass2(Nf,Kn(2)); 

mesh = [25,50,100,200,400]';
ratio =2*[16,8,4,2,1]';
Error = zeros(1,2);


for i = 1:5
    [x, y11, y12, y13, y14] = ass2(mesh(i),Kn(1));
    [x, y21, y22, y23, y24] = ass2(mesh(i),Kn(2));
    k = ratio(i);
    Error(i,1) = max(abs(y11f(1:k:Nf+1)-y11));
    Error(i,2) = max(abs(y21f(1:k:Nf+1)-y21));
end

p1 = polyfit(log(1./mesh),log(Error(:,1)),1);
p2 = polyfit(log(1./mesh),log(Error(:,2)),1);



figure(1);
loglog(mesh,Error(:,1),'--*b','LineWidth',1)
hold on;
loglog(mesh,Error(:,2),'--or','LineWidth',1);
%loglog(mesh,Error(:,3),'-^b','LineWidth',1)
%loglog(mesh,Error(:,4),'-*m','LineWidth',1)
grid on;

hold off;
title('Part 2 (EOC for velocity)','FontSize',13,'Color','k', 'Interpreter', 'latex');
xlabel('$N$','FontSize',13,'Color','k', 'Interpreter', 'latex');
ylabel('$\epsilon_{h}$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
legend('Kn=0.068','Kn=0.5');








