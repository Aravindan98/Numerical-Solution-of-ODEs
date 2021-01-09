  % Stability regions for Adams-Bashforth methods 10
  theta = linspace(0,2*pi,200);
  z = exp(1i*theta); r = z-1;
  s = 1; plot(r./s, 'b')  
  hold on;                                % order 1
  s = (3-1./z)/2; plot(r./s, 'g')                         % order 2
  s = (23-16./z+5./z.^2)/12; plot(r./s, 'r'); % order 3
  s = (55-59./z+37./(z.^2)-9./(z.^3))/24; plot(r./s, 'y'); 
hold on;

xlabel('$Re(z)$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$Im(x)$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
title('Adams-Bashforth','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
axis([-3 1 -2 2]), axis square