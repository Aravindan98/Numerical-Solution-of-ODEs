  % Stability regions for Adams-Moulton methods 3-6
  theta = linspace(0,2*pi,50);
  z = exp(1i*theta); 
  r = (1-z.^(-2));
  s = (1+4*(1./z)+1*(z.^(-2)))/3; 
  plot(r./s, 'r');                    % order 
  %s = (9*z+19-5./z+1./z.^2)/24; plot(r./s, 'c')           % order 4 
  %s = (251*z+646-264./z+106./z.^2-19./z.^3)/720; plot(r./s, 'm')    % 5
  %d = 1-1./z;
  %s = 1-d/2-d.^2/12-d.^3/24-19*d.^4/720-3*d.^5/160; plot(d./s, 'k') % 6
%hold on

xlabel('$Re(z)$','FontSize',13,'Color','k', 'Interpreter', 'latex')
ylabel('$Im(x)$','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
title('Question 2','FontSize',13,'FontWeight','bold','Color','k', 'Interpreter', 'latex')
axis([-7 1 -4 4]), axis square
