% we need \deltai0,omegai0,E'qi0,E'di0,Efdi0,Idi0,Iqi0,Vi0,thetai0,Tmi0 and
% Ybus to calculate the system matrices.
%The initial values calculated
%set the swing bus as bus one, namely bus name11as BusNo1.
% Put bus number 1 as BusNo.3
delta0=[0.7976 0.7610 0.9334 0.5939];
Id0=[5.9958 6.0224 5.7428 6.1174];
Iq0=[3.8886 3.7748 3.8597 3.7171];
Vd0=[0.73716 0.7155 0.7317 0.7045];
Vg0=[0.71937 0.7128 0.7249 0.7237];
Ed0=[0.5016 0.48695 0.4979 0.4795];
Eq0=[0.9183 0.9126 0.9155 0.9266];
Efd0=[1.9196 1.9183 1.8745 1.9482];
Tm0=[7.23122 7.0140 7.0135 7.0142];
%initial value offered in table
V0=[1.03 1.01 1.03 1.01 1.0108 0.9875 1.0095 0.9850 0.9761 0.9716];
theta0=[ 0 -1.5040*pi/180 8.2154*pi/180 -10.2051*pi/180 3.6615*pi/180 -6.2433*pi/180 -4.6977*pi/180 -14.9443*pi/180 -14.4194*pi/180 -23.2922*pi/180];
PG=[7.2172 7 7 7 0 0 0 0 0 0];
QG=[1.4466 1.5920 1.3386 1.8083 0 0 0 0 0 0];
PL=[0 0 0 0 0 0 0 0 11.59 15.75];
GL=[0 0 0 0 0 0 0 0 -0.7350 -0.8990];
%machine data in Table 8.9
X1=[0.022 0.022 0.022 0.022];
Rs=[0.00028 0.00028 0.00028 0.00028];
XXd=[0.2 0.2 0.2 0.2];
Xd=[0.033 0.033 0.033 0.033];
Td0=[8 8 8 8];
XXq=[0.19 0.19 0.19 0.19];
Xq=[0.061 0.061 0.061 0.061];
Tq0=[0.4 0.4 0.4 0.4];
H=[63 54 54 63];
D=[0 0 0 0 ];
omegas=2*pi*60;
M=2*H/omegas;
%Excitation system data in Table8.10
KA=[200 200 200 200];
TA=[0.1 0.1 0.1 0.1];
%Ybus matrix
%set the swing bus as Bus 1
Yb=zeros(10,10);
Yb(3,5)=-(1/(0.001+0.012i));Yb(5,3)=Yb(3,5);
Yb(2,6)=-(1/(0.001+0.012i));Yb(6,2)=Yb(2,6);
Yb(9,10)=-(1/((0.022/3)+(0.22i/3)));Yb(10,9)=Yb(9,10);
Yb(9,6)=-(1/((0.002/2)+(0.02i/2))); Yb(6,9)=Yb(9,6);
Yb(1,7)=-(1/(0.001+0.012i));Yb(7,1)=Yb(1,7);
Yb(4,8)=-(1/(0.001+0.012i));Yb(8,4)=Yb(4,8);
Yb(10,8)=-(1/((0.002/2)+(0.02i/2))); Yb(8,10)=Yb(10,8);
Yb(5,6)=-(1/((0.005/2)+(0.05i/2))); Yb(6,5)=Yb(5,6);
Yb(7,8)=-(1/((0.005/2)+(0.05i/2)));Yb(8,7)=Yb(7,8);
Yb(1,1)=-Yb(1,7);
Yb(2,2)=-Yb(2,6);
Yb(3,3)=-Yb(3,5);
Yb(4,4)=-Yb(4,8);
Yb(5,5)=-(Yb(3,5)+Yb(5,6));
Yb(6,6)=-(Yb(2,6)+Yb(9,6)+Yb(5,6));
Yb(7,7)=-(Yb(1,7)+Yb(7,8));
Yb(8,8)=-(Yb(4,8)+Yb(10,8)+Yb(7,8));
Yb(9,9)=-(Yb(9,10)+Yb(9,6));
Yb(10,10)=-(Yb(9,10)+Yb(10,8));
%get the magnitude and angle of Yb.
Y=zeros(10,10);
alpha=zeros(10,10);
for i=1:10
for j=1:10
Y(i,j)=abs(Yb(i,j));
alpha(i,j)=angle(Yb(i,j));
end
end
%dEfdi=-Efdi/TAi+(KAi/TAi)(Vrefi-Vi)as a first order exciter model
% namely DeltadEfdi=-DeltaEfdi/TAi+(KAi/TAi)(DeltaVrefi-DeltaVi)
%calculate A1,B1,B2,E1,C1,D1,D2
A1=zeros(5,5,4);
B1=zeros(5,2,4);
B2=zeros(5,2,4);
E1=zeros(5,2,4);
C1=zeros(2,5,4);
D1=zeros(2,2,4);
D2=zeros(2,2,4);
for i=1:4;
A1(:,:,i)=[0 1 0 0 0;
0 -D(i)/M(i) -Iq0(i)/M(i) -Id0(i)/M(i) 0;
0 0 -1/Td0(i) 0 1/Td0(i);
0 0 0 -1/Tq0(i) 0;
0 0 0 0 -1/TA(i);];
B1(:,:,i)=[0 0;
(Iq0(i)*(Xd(i)-Xq(i))-Ed0(i))/M(i) (Id0(i)*(Xd(i)-Xq(i))-Eq0(i))/M(i);
-(XXd(i)-Xd(i))/Td0(i) 0;
0 (XXq(i)-Xq(i))/Tq0(i);
0 0];
B2(:,:,i)=[0 0;
0 0;
0 0;
0 0;
0 -KA(i)/TA(i)];
E1(:,:,i)=[0 0;
1/M(i) 0;
0 0;
0 0;
0 KA(i)/TA(i)];
C1(:,:,i)=[-V0(i)*cos(delta0(i)-theta0(i)) 0 0 1 0;
V0(i)*sin(delta0(i)-theta0(i)) 0 1 0 0];
D1(:,:,i)=[-Rs(i) Xq(i)
-Xd(i) -Rs(i)];
D2(:,:,i)=[V0(i)*cos(delta0(i)-theta0(i)) -sin(delta0(i)-theta0(i));
-V0(i)*sin(delta0(i)-theta0(i)) -cos(delta0(i)-theta0(i))];
end
%make these block diagonal matrices
AA1 = blkdiag(A1(:,:,1),A1(:,:,2),A1(:,:,3),A1(:,:,4));
BB1 = blkdiag(B1(:,:,1),B1(:,:,2),B1(:,:,3),B1(:,:,4));
BB2 = blkdiag(B2(:,:,1),B2(:,:,2),B2(:,:,3),B2(:,:,4));
EE1 = blkdiag(E1(:,:,1),E1(:,:,2),E1(:,:,3),E1(:,:,4));
CC1 = blkdiag(C1(:,:,1),C1(:,:,2),C1(:,:,3),C1(:,:,4));
DD1 = blkdiag(D1(:,:,1),D1(:,:,2),D1(:,:,3),D1(:,:,4));
DD2 = blkdiag(D2(:,:,1),D2(:,:,2),D2(:,:,3),D2(:,:,4));
%then we try to get the network equation aboutC2,D3,D4,D5,D6 and D7
C2=zeros(2,5,4);
D3=zeros(2,2,4);
for i=1:4;
C2(:,:,i)=[Id0(i)*V0(i)*cos(delta0(i)-theta0(i))-Iq0(i)*V0(i)*sin(delta0(i)-theta0(i)) 0 0 0 0;
-Id0(i)*V0(i)*sin(delta0(i)-theta0(i))-Iq0(i)*V0(i)*cos(delta0(i)-theta0(i)) 0 0 0 0];
D3(:,:,i)=[V0(i)*sin(delta0(i)-theta0(i)) V0(i)*cos(delta0(i)-theta0(i));
V0(i)*cos(delta0(i)-theta0(i)) -V0(i)*sin(delta0(i)-theta0(i))];
end
%make these block diagonal matrices
CC2=blkdiag(C2(:,:,1),C2(:,:,2),C2(:,:,3),C2(:,:,4));
DD3=blkdiag(D3(:,:,1),D3(:,:,2),D3(:,:,3),D3(:,:,4));
%calculate D4,D5
a=zeros(1,4);
b=zeros(1,4);
c=zeros(1,4);
d=zeros(1,4);
x1=zeros(1,4);
x2=zeros(1,4);
x3=zeros(1,4);
x4=zeros(1,4);
x5=zeros(1,4);
x6=zeros(1,4);
x7=zeros(1,4);
x8=zeros(1,4);
for i=1:4;
x1(i)=Id0(i)*sin(delta0(i)-theta0(i));
x2(i)=-Id0(i)*V0(i)*cos(delta0(i)-theta0(i));
x3(i)=Iq0(i)*cos(delta0(i)-theta0(i));
x4(i)=Iq0(i)*V0(i)*sin(delta0(i)-theta0(i));
x5(i)=Id0(i)*cos(delta0(i)-theta0(i));
x6(i)=Id0(i)*V0(i)*sin(delta0(i)-theta0(i));
x7(i)=-Iq0(i)*sin(delta0(i)-theta0(i));
x8(i)=Iq0(i)*V0(i)*cos(delta0(i)-theta0(i));
for k=1:10;
if k==i;
a(i)=a(i)-V0(k)*Y(i,k)*cos(theta0(i)-theta0(k)-alpha(i,k));
b(i)=b(i)+0;
c(i)=c(i)-V0(k)*Y(i,k)*sin(theta0(i)-theta0(k)-alpha(i,k));
d(i)=d(i)-0;
end
if k~=i;
a(i)=a(i)-V0(k)*Y(i,k)*cos(theta0(i)-theta0(k)-alpha(i,k));
b(i)=b(i)+V0(i)*V0(k)*Y(i,k)*sin(theta0(i)-theta0(k)-alpha(i,k));
c(i)=c(i)-V0(k)*Y(i,k)*sin(theta0(i)-theta0(k)-alpha(i,k));
d(i)=d(i)-V0(i)*V0(k)*Y(i,k)*cos(theta0(i)-theta0(k)-alpha(i,k));
end
end
end
D41=[x2(1)+x4(1)+b(1) x1(1)+x3(1)+a(1)-V0(1)*Y(1,1)*cos(theta0(1)-theta0(1)-alpha(1,1)) 0 0 0 0 0 0;
x6(1)+x8(1)+d(1) x5(1)+x7(1)+c(1)-V0(1)*Y(1,1)*sin(theta0(1)-theta0(1)-alpha(1,1)) 0 0 0 0 0 0];
D42=[0 0 x2(2)+x4(2)+b(2) x1(2)+x3(2)+a(2)-V0(2)*Y(2,2)*cos(theta0(2)-theta0(2)-alpha(2,2)) 0 0 0 0;
0 0 x6(2)+x8(2)+d(2) x5(2)+x7(2)+c(2)-V0(2)*Y(2,2)*sin(theta0(2)-theta0(2)-alpha(2,2)) 0 0 0 0];
D43=[0 0 0 0 x2(3)+x4(3)+b(3) x1(3)+x3(3)+a(3)-V0(3)*Y(3,3)*cos(theta0(3)-theta0(3)-alpha(3,3)) 0 0;
0 0 0 0 x6(3)+x8(3)+d(3) x5(3)+x7(3)+c(3)-V0(3)*Y(3,3)*sin(theta0(3)-theta0(3)-alpha(3,3)) 0 0];
D44=[0 0 0 0 0 0 x2(4)+x4(4)+b(4) x1(4)+x3(4)+a(4)-V0(4)*Y(4,4)*cos(theta0(4)-theta0(4)-alpha(4,4));
0 0 0 0 0 0 x6(4)+x8(4)+d(4) x5(4)+x7(4)+c(4)-V0(4)*Y(4,4)*sin(theta0(4)-theta0(4)-alpha(4,4))];
D51=[0 0 0 0 -V0(1)*V0(7)*Y(1,7)*sin(theta0(1)-theta0(7)-alpha(1,7)) -V0(1)*Y(1,7)*cos(theta0(1)-theta0(7)-alpha(1,7)) 0 0 0 0 0 0;
0 0 0 0 V0(1)*V0(7)*Y(1,7)*cos(theta0(1)-theta0(7)-alpha(1,7)) -V0(1)*Y(1,7)*sin(theta0(1)-theta0(7)-alpha(1,7)) 0 0 0 0 0 0];%1,7
D52=[0 0 -V0(2)*V0(6)*Y(2,6)*sin(theta0(2)-theta0(6)-alpha(2,6)) -V0(2)*Y(2,6)*cos(theta0(2)-theta0(6)-alpha(2,6)) 0 0 0 0 0 0 0 0;
0 0 V0(2)*V0(6)*Y(2,6)*cos(theta0(2)-theta0(6)-alpha(2,6)) -V0(2)*Y(2,6)*sin(theta0(2)-theta0(6)-alpha(2,6)) 0 0 0 0 0 0 0 0];%2,6
D53=[-V0(3)*V0(5)*Y(3,5)*sin(theta0(3)-theta0(5)-alpha(3,5)) -V0(3)*Y(3,5)*cos(theta0(3)-theta0(5)-alpha(3,5)) 0 0 0 0 0 0 0 0 0 0;
V0(3)*V0(5)*Y(3,5)*cos(theta0(3)-theta0(5)-alpha(3,5)) -V0(3)*Y(3,5)*sin(theta0(3)-theta0(5)-alpha(3,5)) 0 0 0 0 0 0 0 0 0 0];%3,5
D54=[0 0 0 0 0 0 -V0(4)*V0(8)*Y(4,8)*sin(theta0(4)-theta0(8)-alpha(4,8)) -V0(4)*Y(4,8)*cos(theta0(4)-theta0(8)-alpha(4,8)) 0 0 0 0;
0 0 0 0 0 0 V0(4)*V0(8)*Y(4,8)*cos(theta0(4)-theta0(8)-alpha(4,8)) -V0(4)*Y(4,8)*sin(theta0(4)-theta0(8)-alpha(4,8)) 0 0 0 0];%4,8
D4=[D41;D42;D43;D44];
D5=[D51;D52;D53;D54];
%calculate D6,D7
e=zeros(1,6);
f=zeros(1,6);
g=zeros(1,6);
h=zeros(1,6);
for i=5:10;
for k=1:10;
if k==i;
e(i-4)=e(i-4)-V0(k)*Y(i,k)*cos(theta0(i)-theta0(k)-alpha(i,k));
f(i-4)=f(i-4)+0;
g(i-4)=g(i-4)-V0(k)*Y(i,k)*sin(theta0(i)-theta0(k)-alpha(i,k));
h(i-4)=h(i-4)-0;
end
if k~=i;
e(i-4)=e(i-4)-V0(k)*Y(i,k)*cos(theta0(i)-theta0(k)-alpha(i,k));
f(i-4)=f(i-4)+V0(i)*V0(k)*Y(i,k)*sin(theta0(i)-theta0(k)-alpha(i,k));
g(i-4)=g(i-4)-V0(k)*Y(i,k)*sin(theta0(i)-theta0(k)-alpha(i,k));
h(i-4)=h(i-4)-V0(i)*V0(k)*Y(i,k)*cos(theta0(i)-theta0(k)-alpha(i,k));
end
end
end
D65=[0 0 0 0 -V0(5)*V0(3)*Y(5,3)*sin(theta0(5)-theta0(3)-alpha(5,3)) -V0(5)*Y(5,3)*cos(theta0(5)-theta0(3)-alpha(5,3)) 0 0;
0 0 0 0 V0(5)*V0(3)*Y(5,3)*cos(theta0(5)-theta0(3)-alpha(5,3)) -V0(5)*Y(5,3)*sin(theta0(5)-theta0(3)-alpha(5,3)) 0 0];%5,3
D66=[0 0 -V0(6)*V0(2)*Y(6,2)*sin(theta0(6)-theta0(2)-alpha(6,2)) -V0(6)*Y(6,2)*cos(theta0(6)-theta0(2)-alpha(6,2)) 0 0 0 0;
0 0 V0(6)*V0(2)*Y(6,2)*cos(theta0(6)-theta0(2)-alpha(6,2)) -V0(6)*Y(6,2)*sin(theta0(6)-theta0(2)-alpha(6,2)) 0 0 0 0];%6,2
D67=[-V0(7)*V0(1)*Y(7,1)*sin(theta0(7)-theta0(1)-alpha(7,1)) -V0(7)*Y(7,1)*cos(theta0(7)-theta0(1)-alpha(7,1)) 0 0 0 0 0 0;
V0(7)*V0(1)*Y(7,1)*cos(theta0(7)-theta0(1)-alpha(7,1)) -V0(7)*Y(7,1)*sin(theta0(7)-theta0(1)-alpha(7,1)) 0 0 0 0 0 0];%7,1
D68=[0 0 0 0 0 0 -V0(8)*V0(4)*Y(8,4)*sin(theta0(8)-theta0(4)-alpha(8,4)) -V0(8)*Y(8,4)*cos(theta0(8)-theta0(4)-alpha(8,4));
0 0 0 0 0 0 V0(8)*V0(4)*Y(8,4)*cos(theta0(8)-theta0(4)-alpha(8,4)) -V0(8)*Y(8,4)*sin(theta0(8)-theta0(4)-alpha(8,4))];%8,4
D69=[0 0 0 0 0 0 0 0;%none
0 0 0 0 0 0 0 0];
D610=[0 0 0 0 0 0 0 0;%none
0 0 0 0 0 0 0 0];
D75=[f(1) e(1)-V0(5)*Y(5,5)*cos(theta0(5)-theta0(5)-alpha(5,5)) -V0(5)*V0(6)*Y(5,6)*sin(theta0(5)-theta0(6)-alpha(5,6)) -V0(5)*Y(5,6)*cos(theta0(5)-theta0(6)-alpha(5,6)) 0 0 0 0 0 0 0 0;
h(1) g(1)-V0(5)*Y(5,5)*sin(theta0(5)-theta0(5)-alpha(5,5)) V0(5)*V0(6)*Y(5,6)*cos(theta0(5)-theta0(6)-alpha(5,6)) -V0(5)*Y(5,6)*sin(theta0(5)-theta0(6)-alpha(5,6)) 0 0 0 0 0 0 0 0];%5,5-5,6
D76=[-V0(6)*V0(5)*Y(6,5)*sin(theta0(6)-theta0(5)-alpha(6,5)) -V0(6)*Y(6,5)*cos(theta0(6)-theta0(5)-alpha(6,5)) f(2) e(2)-V0(6)*Y(6,6)*cos(theta0(6)-theta0(6)-alpha(6,6)) 0 0 0 0 -V0(6)*V0(9)*Y(6,9)*sin(theta0(6)-theta0(9)-alpha(6,9)) -V0(5)*Y(6,9)*cos(theta0(6)-theta0(9)-alpha(6,9)) 0 0;
V0(6)*V0(5)*Y(6,5)*cos(theta0(6)-theta0(5)-alpha(6,5)) -V0(6)*Y(6,5)*sin(theta0(6)-theta0(5)-alpha(6,5)) h(2) g(2)-V0(6)*Y(6,6)*sin(theta0(6)-theta0(6)-alpha(6,6)) 0 0 0 0 V0(6)*V0(9)*Y(6,9)*cos(theta0(6)-theta0(9)-alpha(6,9)) -V0(5)*Y(6,9)*sin(theta0(6)-theta0(6)-alpha(6,9)) 0 0];%6,5-6,6-6,9
D77=[0 0 0 0 f(3) e(3)-V0(7)*Y(7,7)*cos(theta0(7)-theta0(7)-alpha(7,7)) -V0(7)*V0(8)*Y(7,8)*sin(theta0(7)-theta0(8)-alpha(7,8)) -V0(7)*Y(7,8)*cos(theta0(7)-theta0(8)-alpha(7,8)) 0 0 0 0;
0 0 0 0 h(3) g(3)-V0(7)*Y(7,7)*sin(theta0(7)-theta0(7)-alpha(7,7)) V0(7)*V0(8)*Y(7,8)*cos(theta0(7)-theta0(8)-alpha(7,8)) -V0(7)*Y(7,8)*sin(theta0(7)-theta0(8)-alpha(7,8)) 0 0 0 0];%7,7-7,8
D78=[0 0 0 0 -V0(8)*V0(7)*Y(8,7)*sin(theta0(8)-theta0(7)-alpha(8,7)) -V0(8)*Y(8,7)*cos(theta0(8)-theta0(7)-alpha(8,7)) f(4) e(4)-V0(8)*Y(8,8)*cos(theta0(8)-theta0(8)-alpha(8,8)) 0 0 -V0(8)*V0(10)*Y(8,10)*sin(theta0(8)-theta0(10)-alpha(8,10)) -V0(8)*Y(8,10)*cos(theta0(8)-theta0(10)-alpha(8,10));
0 0 0 0 V0(8)*V0(7)*Y(8,7)*cos(theta0(8)-theta0(7)-alpha(8,7)) -V0(8)*Y(8,7)*sin(theta0(8)-theta0(7)-alpha(8,7)) h(4) g(4)-V0(8)*Y(8,8)*sin(theta0(8)-theta0(8)-alpha(8,8)) 0 0 V0(8)*V0(10)*Y(8,10)*cos(theta0(8)-theta0(10)-alpha(8,10)) -V0(8)*Y(8,10)*sin(theta0(8)-theta0(10)-alpha(8,10))];%8,7-8,8-8,10
D79=[0 0 -V0(9)*V0(6)*Y(9,6)*sin(theta0(9)-theta0(6)-alpha(9,6)) -V0(9)*Y(9,6)*cos(theta0(9)-theta0(6)-alpha(9,6)) 0 0 0 0 f(5) e(5)-V0(9)*Y(9,9)*cos(theta0(9)-theta0(9)-alpha(9,9)) -V0(9)*V0(10)*Y(9,10)*sin(theta0(9)-theta0(10)-alpha(9,10)) -V0(9)*Y(9,10)*cos(theta0(9)-theta0(10)-alpha(9,10));
0 0 V0(9)*V0(6)*Y(9,6)*cos(theta0(9)-theta0(6)-alpha(9,6)) -V0(9)*Y(9,6)*sin(theta0(9)-theta0(6)-alpha(9,6)) 0 0 0 0 h(5) g(5)-V0(9)*Y(9,9)*sin(theta0(9)-theta0(9)-alpha(9,9)) V0(9)*V0(10)*Y(9,10)*cos(theta0(9)-theta0(10)-alpha(9,10)) -V0(9)*Y(9,10)*sin(theta0(9)-theta0(10)-alpha(9,10))];%9,6-9,9-9,10
D710=[0 0 0 0 0 0 -V0(10)*V0(8)*Y(10,8)*sin(theta0(10)-theta0(8)-alpha(10,8)) -V0(10)*Y(10,8)*cos(theta0(10)-theta0(8)-alpha(10,8)) -V0(10)*V0(9)*Y(10,9)*sin(theta0(10)-theta0(9)-alpha(10,9)) -V0(10)*Y(10,9)*cos(theta0(10)-theta0(9)-alpha(10,9)) f(6) e(6)-V0(10)*Y(10,10)*cos(theta0(10)-theta0(10)-alpha(10,10));
0 0 0 0 0 0 V0(10)*V0(8)*Y(10,8)*cos(theta0(10)-theta0(8)-alpha(10,8)) -V0(10)*Y(10,8)*sin(theta0(10)-theta0(8)-alpha(10,8)) V0(10)*V0(9)*Y(10,9)*cos(theta0(10)-theta0(9)-alpha(10,9)) -V0(10)*Y(10,9)*sin(theta0(10)-theta0(9)-alpha(10,9)) h(6) g(6)-V0(10)*Y(10,10)*sin(theta0(10)-theta0(10)-alpha(10,10))];%10,8-10,9-10,10
D6=[D65;D66;D67;D68;D69;D610];
D7=[D75;D76;D77;D78;D79;D710];
%simplify equations
Anew=AA1-BB1/DD1*CC1;
Bnew=BB2-BB1/DD1*DD2;
XX=zeros(20,12);
Bnn=[Bnew XX];
K1=D4-DD3/DD1*DD2;
K2=CC2-DD3/DD1*CC1;
Cnew=K2;
X2=zeros(12,20);
Cnn=[Cnew;
X2];
Dnew=[K1 D5;
D6 D7];
Asys=Anew-Bnn*inv(Dnew)*Cnn;
%get the right eigenvectors, left eigenvectors, and eigenvalues of Asys
[V,D] = eig(Asys);
W=inv(V)';
eigenvalues=zeros(1,20);
for i=1:20
eigenvalues(i)=D(i,i);
end

eigenvalues

%get real part and imaginary part of the eigenvalues
sigma=zeros(1,20);
omega=zeros(1,20);
for i=1:20;
sigma(i)=real(eigenvalues(i));
omega(i)=abs(imag(eigenvalues(i)));
end
%get the frequency and damping ratio
frequency=zeros(1,20);
dampingratio=zeros(1,20);
for i=1:20
frequency(i)=omega(i)/(2*pi);
dampingratio(i)=(-1)*sigma(i)*((sigma(i)^2+omega(i)^2)^(-0.5));
end
%get the participation factor
WV=zeros(1,20);
for i=1:20;
WV(i)=W(:,i)'*V(:,i);
end
paticipationfactor=zeros(20,20);
for i=1:20;
for j=1:20;
paticipationfactor(i,j)=abs(W(j,i)*V(j,i)/WV(i));
end
end
paticipationfactor=paticipationfactor';
%max of row
maxd=zeros(1,20);
for i=1:20;
maxd(i)=max(paticipationfactor(i,:));
end
%normalize to 1
for i=1:20;
for j=1:20;
paticipationfactor(i,j)=1/maxd(i)*paticipationfactor(i,j);
end
end
% to calculate Rhi we need to get the b and c first
Bi=EE1;
Bb=zeros(20,20);
for i=1:20
for j=1:8;
Bb(i,j)=EE1(i,j);
end
end
Cc=zeros(20,20);
Ci=zeros(4,20);
Ci(1,2)=1;
Ci(2,7)=1;
Ci(3,12)=1;
Ci(4,17)=1;
for i=1:4;
for j=1:20;
Cc(i,j)=Ci(i,j);
end
end
Rh=W'*Bb*Cc'*V;
%calculate phi,a and tau
phi=zeros(20,20);
ai=zeros(20,20);
tau=zeros(20,20);
for i=1:20;
for j=1:20;
phi(i,j)=(pi-angle(Rh(i,j)))/2;
ai(i,j)=(1+sin(phi(i,j)))/(1-sin(phi(i,j)));
tau(i,j)=1/(omega(i)*(ai(i,j)^0.5));
end
end

xlswrite('PSS_design_yijun.xlsx',ai,1) 
xlswrite('PSS_design_yijun.xlsx',tau,2)

% modeshape910 = [V(1,9) V(6,9) V(11,9) V(16,9)];
%modeshape1112 = [V(1,11) V(6,11) V(11,11) V(16,11)];
modeshape1314 = [V(1,13) V(6,13) V(11,13) V(16,13)];
figure
compass(modeshape1314)