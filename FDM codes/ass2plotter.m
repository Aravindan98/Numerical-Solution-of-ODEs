close all;
clear all;


N=200;
k=[0.068,0.1,0.5];
[x, y11, y21, y31, y41] = ass2(N,k(1));
[x, y12, y22, y32, y42] = ass2(N,k(2));
[x, y13, y23, y33, y43] = ass2(N,k(3));

figure(1);
grid on;
hold on;
plot(x,y11,'-r','LineWidth',2);
plot(x,y12,'-g','LineWidth',2);
plot(x,y13,'-b','LineWidth',2);
hold off;
title('Part 1(a)','FontSize',13,'Color','k', 'Interpreter', 'latex')
xlabel('$y$','FontSize',13,'Color','k', 'Interpreter', 'latex');
ylabel('$v_{x}$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
legend('Kn=0.068','Kn=0.1','Kn=0.5');

figure(2);
grid on;
hold on;
plot(x,y21,'-r','LineWidth',2);
plot(x,y22,'-g','LineWidth',2);
plot(x,y23,'-b','LineWidth',2);
hold off;
title('Part 1(b)','FontSize',13,'Color','k', 'Interpreter', 'latex')
xlabel('$y$','FontSize',13,'Color','k', 'Interpreter', 'latex');
ylabel('$\sigma$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
legend('Kn=0.068','Kn=0.1','Kn=0.5');

figure(3);
grid on;
hold on;
plot(x,y31,'-r','LineWidth',2);
plot(x,y32,'-g','LineWidth',2);
plot(x,y33,'-b','LineWidth',2);
hold off;
title('Part 1(c)','FontSize',13,'Color','k', 'Interpreter', 'latex')
xlabel('$y$','FontSize',13,'Color','k', 'Interpreter', 'latex');
ylabel('$\theta$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
legend('Kn=0.068','Kn=0.1','Kn=0.5');

figure(4);
grid on;
hold on;
plot(x,y41,'-r','LineWidth',2);
plot(x,y42,'-g','LineWidth',2);
plot(x,y43,'-b','LineWidth',2);
hold off;
title('Part 1(d)','FontSize',13,'Color','k', 'Interpreter', 'latex');
xlabel('$y$','FontSize',13,'Color','k', 'Interpreter', 'latex');
ylabel('$q_{y}$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex');
legend('Kn=0.068','Kn=0.1','Kn=0.5');


