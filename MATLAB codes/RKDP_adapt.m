function [y] = RKDP_adapt(fcn, tspan, y0, n)
% classical fourth order method
% Solve the initial value problem
% y’ = fcn(t,y), t0 <= t <= b, y(t0)=y0

tn_1 = tspan(1);
%t_end = tspan(2);
h = (tspan(2)-tspan(1))/n; % step size.
%t = t0:h:t_end;
tol = 1e-3;
%y = zeros(length(y0),length(t)); % solution vector initialization
%ys = zeros(length(y0),length(t)); % solution vector initialization
yn_1 = y0;
ysn_1=y0;
hn_1 = h;
%y(:,1) = y0; % initial condition
%ys(:,1) = y0; % initial condition
% J. R. Dormand and P. J. Prince coeffcients
y = y0;
t = tspan(1);
A = [0, 0, 0, 0, 0, 0, 0;
    1/5, 0, 0, 0, 0, 0, 0;
    3/40, 9/40, 0, 0, 0, 0, 0;
    44/45, -56/15, 32/9, 0, 0, 0, 0;
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0;
    9017/3168, -355/33, 46732/5247, 49/176 , -5103/18656, 0, 0;
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
    

b1 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]';
b2 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]';

c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]';

%for i = 2:length(t

while tn_1<tspan(2)
	tn = tn_1;
	yn = yn_1;
	ysn = ysn_1;
	h = hn_1;
	disp(h);
	disp(tn_1);
    k_1 = feval(fcn,tn,yn);
    k_2 = feval(fcn,tn+c(2)*h,yn+h*A(2,1)*k_1);
    k_3 = feval(fcn,tn+c(3)*h,yn+h*A(3,1)*k_1+h*A(3,2)*k_2);
    k_4 = feval(fcn,tn+c(4)*h,yn+h*A(4,1)*k_1+h*A(4,2)*k_2+h*A(4,3)*k_3);
    k_5 = feval(fcn,tn+c(5)*h,yn+h*A(5,1)*k_1+h*A(5,2)*k_2+h*A(5,3)*k_3+h*A(5,4)*k_4);
    k_6 = feval(fcn,tn+c(6)*h,yn+h*A(6,1)*k_1+h*A(6,2)*k_2+h*A(6,3)*k_3+h*A(6,4)*k_4+h*A(6,5)*k_5);
    k_7 = feval(fcn,tn+c(7)*h,yn+h*A(7,1)*k_1+h*A(7,2)*k_2+h*A(7,3)*k_3+h*A(7,4)*k_4+h*A(7,5)*k_5+h*A(7,6)*k_6);

    %y(:,i) = yn+h*b1(1)*k_1+h*b1(2)*k_2+h*b1(3)*k_3+h*b1(4)*k_4+h*b1(5)*k_5+h*b1(6)*k_6+h*b1(7)*k_7;
    %ys(:,i) = ys(:,i-1)+h*b2(1)*k_1+h*b2(2)*k_2+h*b2(3)*k_3+h*b2(4)*k_4+h*b2(5)*k_5+h*b2(6)*k_6+h*b2(7)*k_7;
    error = max(abs(yn+h*b1(1)*k_1+h*b1(2)*k_2+h*b1(3)*k_3+h*b1(4)*k_4+h*b1(5)*k_5+h*b1(6)*k_6+h*b1(7)*k_7 - (ysn+h*b2(1)*k_1+h*b2(2)*k_2+h*b2(3)*k_3+h*b2(4)*k_4+h*b2(5)*k_5+h*b2(6)*k_6+h*b2(7)*k_7)));
 
    if error < tol
        hn_1 = max(min(0.9*h*((tol/error)^0.25),0.02),1e-3);
    else
        hn_1 = min(max(0.9*h*((tol/error)^0.2),1e-3),0.02);
    end
    
	k_1 = feval(fcn,tn,yn);
    k_2 = feval(fcn,tn+c(2)*hn_1,yn+hn_1*A(2,1)*k_1);
    k_3 = feval(fcn,tn+c(3)*hn_1,yn+hn_1*A(3,1)*k_1+hn_1*A(3,2)*k_2);
    k_4 = feval(fcn,tn+c(4)*hn_1,yn+hn_1*A(4,1)*k_1+hn_1*A(4,2)*k_2+hn_1*A(4,3)*k_3);
    k_5 = feval(fcn,tn+c(5)*hn_1,yn+hn_1*A(5,1)*k_1+hn_1*A(5,2)*k_2+hn_1*A(5,3)*k_3+hn_1*A(5,4)*k_4);
    k_6 = feval(fcn,tn+c(6)*hn_1,yn+hn_1*A(6,1)*k_1+hn_1*A(6,2)*k_2+hn_1*A(6,3)*k_3+hn_1*A(6,4)*k_4+hn_1*A(6,5)*k_5);
    k_7 = feval(fcn,tn+c(7)*hn_1,yn+hn_1*A(7,1)*k_1+hn_1*A(7,2)*k_2+hn_1*A(7,3)*k_3+hn_1*A(7,4)*k_4+hn_1*A(7,5)*k_5+hn_1*A(7,6)*k_6);
	
	yn_1=yn+hn_1*b1(1)*k_1+hn_1*b1(2)*k_2+hn_1*b1(3)*k_3+hn_1*b1(4)*k_4+hn_1*b1(5)*k_5+hn_1*b1(6)*k_6+hn_1*b1(7)*k_7;
	ysn_1 = ysn+hn_1*b2(1)*k_1+hn_1*b2(2)*k_2+hn_1*b2(3)*k_3+hn_1*b2(4)*k_4+hn_1*b2(5)*k_5+hn_1*b2(6)*k_6+hn_1*b2(7)*k_7;
	tn_1 = tn + hn_1;
	y = [y,yn_1];
	%t = [t,tn_1];
end
%error = abs(y-ys);
end