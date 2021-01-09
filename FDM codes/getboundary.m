function [Ba,Bb,alpha]=util(g,y,N,Ba,Bb,mm)
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