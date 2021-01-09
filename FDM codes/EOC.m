% ASR NODES 2020
% The non-linear mid-point method
% SOLVES y'(x)=q(x,y) 
% g(y(a),g(b)) = 0 
%Ex 6
%EOC
format long

Nf = 800;
[xf, y1f, y2f, y3f, y4f] = FDEx4(Nf, 10);



mesh = [25,50,100,200]';
ratio =2*[16,8,4,2]';
Error = zeros(4,4);


for i = 1:4
    [x, y1, y2, y3, y4] = FDEx4(mesh(i), 10);
    k = ratio(i);
    Error(i,1) = max(abs(y1f(1:k:Nf+1)-y1)); % l_inf error in y1
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
hold on
loglog(mesh,Error(:,2),'-sr','LineWidth',1)
loglog(mesh,Error(:,3),'-^b','LineWidth',1)
loglog(mesh,Error(:,4),'-*m','LineWidth',1)
grid on;

xlabel('$N$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$\epsilon_{h}$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')


