function [t,y] = classicalRK_4(fcn, tspan, y0, n)
% classical fourth order method
% Solve the initial value problem
% y’ = fcn(t,y), t0 <= t <= b, y(t0)=y0

t0 = tspan(1);
t_end = tspan(2);
h = (t_end-t0)/n; % step size.
t = t0:h:t_end;
y = zeros(length(y0),length(t)); % solution vector initialization

y(:,1) = y0; % initial condition

% example
 for i = 2:length(t)
    k_1 = feval(fcn,t(i-1),y(:,i-1));
    k_2 = feval(fcn,t(i-1)+0.5*h,y(:,i-1)+0.5*h*k_1);
    k_3 = feval(fcn,t(i-1)+0.5*h,y(:,i-1)+0.5*h*k_2);
    k_4 = feval(fcn,t(i-1)+h,y(:,i-1)+h*k_3);
    
    y(:,i) = y(:,i-1) + (1/6)*h*(k_1+2*k_2+2*k_3+k_4);
end

end
