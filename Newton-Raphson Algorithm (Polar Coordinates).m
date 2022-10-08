clc
clear
close all
%% 计算节点导纳矩阵
G(1,1)=10.2;B(1,1)=-31.5;G(1,2)=-1.2;B(1,2)=4.0;G(1,3)=-1.5;
B(1,3)=5.0;G(1,4)=-2.5;B(1,4)=7.5;G(1,5)=-5.000;B(1,5)=15.000;
 
G(2,1)=-1.2;B(2,1)=4.0; G(2,2)=10.4;B(2,2)=-31.7;G(2,3)=-8.0;
B(2,3)=24.0;G(2,4)=0;B(2,4)=0;G(2,5)=-1.2;B(2,5)=3.7;
 
G(3,1)=-1.5;B(3,1)=5.0;G(3,2)=-8.0;B(3,2)=24.0;G(3,3)=10.7;B(3,3)=-32.7; G(3,4)=-1.2;B(3,4)=3.7;G(3,5)=0;B(3,5)=0;
 
G(4,1)=-2.500;B(4,1)=7.500;G(4,2)=0;B(4,2)=0;G(4,3)=-1.2;B(4,3)=3.7;G(4,4)=3.7;B(4,4)=-11.2;G(4,5)=0;B(4,5)=0;
 
G(5,1)=-5.0;B(5,1)=15.0;G(5,2)=-1.2;B(5,2)=3.7;G(5,3)=0;B(5,3)=0;G(5,4)=0;B(5,4)=0;G(5,5)=6.2;B(5,5)=-18.7;
 
Y=G+j*B;
%% 设置各节点电压初值 
delt(1)=0;delt(2)=0;delt(3)=0;delt(4)=0;
u(1)=1.0;u(2)=1.0;u(3)=1.0;u(4)=1.0;
%% 输入原始数据
p(1)=0.20; q(1)=0.20; p(2)=-0.45; q(2)=-0.15; 
p(3)=-0.40; q(3)=-0.05; p(4)=-0.60; q(4)=-0.10;
d(1,4)=0; d(4,1)=0; d(1,5)=0;d(5,1)=0; 
%% 牛顿-拉夫逊算法(极坐标) 
k=0;precision=1; 
k,delt,u 
N1=4; 
while precision>0.00001 
    delt(5)=0;u(5)=1.05;
for m=1:N1 
for n=1:N1+1 
    pt(n)=u(m)*u(n)*(G(m,n)*cos(delt(m)-delt(n))+B(m,n)*sin(delt(m)-delt(n))); 
    qt(n)=u(m)*u(n)*(G(m,n)*sin(delt(m)-delt(n))-B(m,n)*cos(delt(m)-delt(n))); 
end 
    pp(m)=p(m)-sum(pt);qq(m)=q(m)-sum(qt); 
end 
pp,qq

% 计算雅可比矩阵对角线元素（i＝j时）
for m=1:N1 
for n=1:N1+1 
h0(n)=u(m)*u(n)*(G(m,n)*sin(delt(m)-delt(n))-B(m,n)*cos(delt(m)-delt(n))); 
n0(n)=-u(m)*u(n)*(G(m,n)*cos(delt(m)-delt(n))+B(m,n)*sin(delt(m)-delt(n))); 
j0(n)=-u(m)*u(n)*(G(m,n)*cos(delt(m)-delt(n))+B(m,n)*sin(delt(m)-delt(n))); 
l0(n)=-u(m)*u(n)*(G(m,n)*sin(delt(m)-delt(n))-B(m,n)*cos(delt(m)-delt(n))); 
end 
H(m,m)=sum(h0)-u(m)^2*(G(m,m)*sin(delt(m)-delt(m))-B(m,m)*cos(delt(m)-delt(m))); 
N(m,m)=sum(n0)-2*u(m)^2*G(m,m)+u(m)^2*(G(m,m)*cos(delt(m)-delt(m)))+B(m,m)*sin(delt(m)-delt(m)); 
J(m,m)=sum(j0)+u(m)^2*(G(m,m)*cos(delt(m)-delt(m))+B(m,m)*sin(delt(m)-delt(m))); 
L(m,m)=sum(l0)+2*u(m)^2*B(m,m)+u(m)^2*(G(m,m)*sin(delt(m)-delt(m)) -B(m,m)*cos(delt(m)-delt(m))); 
end
for m=1:N1 
JJ(2*m-1,2*m-1)=H(m,m);JJ(2*m-1,2*m)=N(m,m); 
JJ(2*m,2*m-1)=J(m,m);JJ(2*m,2*m)=L(m,m); 
end 
% 计算雅可比矩阵非对角线元素（i≠j时）
for m=1:N1 
for n=1:N1 
if m==n 
else 
H(m,n)=-u(m)*u(n)*(G(m,n)*sin(delt(m)-delt(n))-B(m,n)*cos(delt(m)-delt(n))); 
J(m,n)=u(m)*u(n)*(G(m,n)*cos(delt(m)-delt(n))+B(m,n)*sin(delt(m)-delt(n))); 
N(m,n)=-J(m,n);L(m,n)=H(m,n); 
end 
end 
end 
for m=1:N1  
for n=1:N1 
if m==n 
else 
JJ(2*m-1,2*n-1)=H(m,n);JJ(2*m-1,2*n)=N(m,n); 
JJ(2*m,2*n-1)=J(m,n);JJ(2*m,2*n)=L(m,n); 
end 
end 
end 
 
for m=1:N1 
PP(2*m-1)=pp(m);PP(2*m)=qq(m); 
end 
uu=-inv(JJ)*PP';precision=max(abs(uu));uu
for n=1:N1 
delt(n)=delt(n)+uu(2*n-1); 
u(n)=u(n)+uu(2*n); 
end 
k=k+1; 
k,delt,u 
end
%% 计算潮流 
for n=1:N1+1 
 
U(n)=u(n)*(cos(delt(n))+j*sin(delt(n))); 
end 
for m=1:N1+1 
I(m)=Y(5,m)*U(m); 
end 
S5=U(5)*sum(conj(I)); 
for n=1:N1+1 
    q4(n)=u(4)*u(n)*(G(4,n)*sin(delt(4)-delt(n))-B(4,n)*cos(delt(4)-delt(n))); 
end 
Q4=sum(q4) 
for m=1:N1+1 
for n=1:N1+1 
S(m,n)=U(m)*(conj(U(m))*conj(d(m,n))+(conj(U(m))-conj(U(n)))*conj(-Y(m,n))); 
end 
end
Y 
JJ 
S 
B 	
pp 
qq 
uu 
U 
k 
Q4 
S5