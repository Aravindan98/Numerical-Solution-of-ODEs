function [x, y1, y2, y3, y4] = FDM(N,k)
mm = 4;
%k = 0.1;
F = 0.23;
h = 1/N;

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
end