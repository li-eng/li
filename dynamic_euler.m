 clear all
tic


% t=0:0.0001:1;%x,y向位移仿真区间
t=0:0.0001:1;%z向位移仿真区间
N=length(t);
% h=0.0001;%x,y向位移步长
h=0.0001;%z向位移步长

%读入z向位移需要的力
 f=xlsread('f_ry_0.0001.xlsx');
%读入x向位移需要的力
% f=xlsread('fy.xlsx');
%读入x向位移需要的力
%  f=xlsread('fx.xlsx');

%分配空间

y1=zeros(1,N);y2=zeros(1,N);y3=zeros(1,N);
dy1=zeros(1,N);dy2=zeros(1,N);dy3=zeros(1,N);
ddy1=zeros(1,N);ddy2=zeros(1,N);ddy3=zeros(1,N);

alpha=zeros(1,N);beta=zeros(1,N);gama=zeros(1,N);
ald=zeros(1,N);bed=zeros(1,N);gad=zeros(1,N);
aldd=zeros(1,N);bedd=zeros(1,N);gadd=zeros(1,N);
Ax=zeros(6,N);
%初值
y1(1)=0;y2(1)=0;y3(1)=0.3825;alpha(1)=0;beta(1)=0.00001;gama(1)=0;
dy1(1)=0;dy2(1)=0;dy3(1)=0;ald(1)=0;bed(1)=0;gad(1)=0;

   
for i=1:N
    f1=f(1,:);f2=f(2,:);f3=f(3,:);f4=f(4,:);f5=f(5,:);f6=f(6,:);
   
    ddy1(i)=ax1(y1(i),y2(i),y3(i),alpha(i),beta(i),gama(i),dy1(i),dy2(i),dy3(i),ald(i),bed(i),gad(i),f1(i),f2(i),f3(i),f4(i),f5(i),f6(i));
    ddy2(i)=ax2(y1(i),y2(i),y3(i),alpha(i),beta(i),gama(i),dy1(i),dy2(i),dy3(i),ald(i),bed(i),gad(i),f1(i),f2(i),f3(i),f4(i),f5(i),f6(i));
    ddy3(i)=ax3(y1(i),y2(i),y3(i),alpha(i),beta(i),gama(i),dy1(i),dy2(i),dy3(i),ald(i),bed(i),gad(i),f1(i),f2(i),f3(i),f4(i),f5(i),f6(i));
    aldd(i)=ax4(y1(i),y2(i),y3(i),alpha(i),beta(i),gama(i),dy1(i),dy2(i),dy3(i),ald(i),bed(i),gad(i),f1(i),f2(i),f3(i),f4(i),f5(i),f6(i));
    bedd(i)=ax5(y1(i),y2(i),y3(i),alpha(i),beta(i),gama(i),dy1(i),dy2(i),dy3(i),ald(i),bed(i),gad(i),f1(i),f2(i),f3(i),f4(i),f5(i),f6(i));
    gadd(i)=ax6(y1(i),y2(i),y3(i),alpha(i),beta(i),gama(i),dy1(i),dy2(i),dy3(i),ald(i),bed(i),gad(i),f1(i),f2(i),f3(i),f4(i),f5(i),f6(i));

    y1(i+1)=y1(i)+dy1(i)*h+0.5*ddy1(i)*h^2;%线性运动迭代
    y2(i+1)=y2(i)+dy2(i)*h+0.5*ddy2(i)*h^2;
    y3(i+1)=y3(i)+dy3(i)*h+0.5*ddy3(i)*h^2;
    alpha(i+1)=alpha(i)+ald(i)*h+0.5*aldd(i)*h^2;
    beta(i+1)=beta(i)+bed(i)*h+0.5*bedd(i)*h^2;
    gama(i+1)=gama(i)+gad(i)*h+0.5*gadd(i)*h^2;
    dy1(i+1)=dy1(i)+ddy1(i)*h;%线运动一阶导迭代
    dy2(i+1)=dy2(i)+ddy2(i)*h;
    dy3(i+1)=dy3(i)+ddy3(i)*h;
    ald(i+1)=ald(i)+aldd(i)*h;%角运动一阶导迭代
    bed(i+1)=bed(i)+bedd(i)*h;
    gad(i+1)=gad(i)+gadd(i)*h;
end
% x=pi/2*t.^2-pi/3*t.^3;
% error=x-theta(1:10001);
plot(t,y3(1:10001))
toc
function Ax=ax(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16)
% 上平台铰链位置
b1=[0.07320508076,-0.07320508076,-0.0365];
b2=[-0.0730508076,-0.07320508076,-0.0365];
b3=[-0.1,-0.02679491924,-0.0365];
b4=[-0.02679491924,0.1,-0.0365];
b5=[0.02679491924,0.1,-0.0365];
b6=[0.1,-0.02679491924,-0.0365];
% 下平台铰链位置
a1=[0.045,-0.16794228634,0.03025];
a2=[-0.045,-0.16794228634,0.03025];
a3=[-0.16794228634,0.045,0.03025];
a4=[-0.12294228634,0.12294228634,0.03025];
a5=[0.12294228634,0.12294228634,0.03025];
a6=[0.16794228634,0.045,0.03025];


%旋转矩阵
sa=sin(alpha);
ca=cos(alpha);
sb=sin(beta);
cb=cos(beta);
sg=sin(gama);
cg=cos(gama);
R=[cb*cg sa*sb*cg-ca*sg ca*sb*cg+sg*sa;cb*sg sa*sb*sg+ca*cg ca*sb*sg-cg*sa;-sb sa*cb ca*cb];
%角运动
tr=R(1,1)+R(2,2)+R(3,3);
mid=(tr-1)/2;
theta=acos(mid);
sx=(R(3,2)-R(2,3))/(2*sin(theta));
sy=(R(1,3)-R(3,1))/(2*sin(theta));
sz=(R(2,1)-R(1,2))/(2*sin(theta));

r1d=-bed*sb*cg-gad*sg*cb;
r2d=ald*ca*sb*sg+bed*sa*cb*sg+gad*sa*sb*sg;
r3d=-ald*sa*cg-gad*sg*ca;
r4d=-ald*sa*cb-bed*sb*ca;
trd=r1d+r2d+r3d+r4d;

omiga=-(3-tr^2+2*tr)^(-0.5)*trd;
omiga1=sx*omiga;
omiga2=sy*omiga;
omiga3=sz*omiga;
w=[omiga1 omiga2 omiga3];

 %直线运动
 p=[y1 y2 y3];Vp=[dy1 dy2 dy3];
 %机构参数
m_dj=0.667417;m_sg=0.072068;m_gz=0.07012;
m11=m_dj+2*m_gz+m_sg;m12=2*m_gz;
c11=0.029811+0.024+0.0365;c12=0.09+0.021;%+0.046
g=[0 0 -9.8];
   %杆的转动惯量
  Ixx=0.003395+0.000391;
% 上平台质量及惯性矩
m=0.8151490458;Ip=[0.00403852037 0 0;0 0.0040385 0;0 0 0.0080758];
AIp=R*Ip*R';
 %静坐标系下的上铰链点位置
 br1=R*b1'; br2=R*b2'; br3=R*b3'; br4=R*b4'; br5=R*b5'; br6=R*b6';%列向量
 % 杆长及其单位向量
  l1=sqrt(p*p'+b1*b1'+a1*a1'-2*p*a1'+2*p*br1-2*br1'*a1');
  l2=sqrt(p*p'+b2*b2'+a2*a2'-2*p*a2'+2*p*br2-2*br2'*a2');
  l3=sqrt(p*p'+b3*b3'+a3*a3'-2*p*a3'+2*p*br3-2*br3'*a3');
  l4=sqrt(p*p'+b4*b4'+a4*a4'-2*p*a4'+2*p*br4-2*br4'*a4');
  l5=sqrt(p*p'+b5*b5'+a5*a5'-2*p*a5'+2*p*br5-2*br5'*a5');
  l6=sqrt(p*p'+b6*b6'+a6*a6'-2*p*a6'+2*p*br6-2*br6'*a6');

  s1=(p+br1'-a1)/l1;%行向量
  s2=(p+br2'-a2)/l2;
  s3=(p+br3'-a3)/l3;
  s4=(p+br4'-a4)/l4;
  s5=(p+br5'-a5)/l5;
  s6=(p+br6'-a6)/l6;
  %单位向量的反对称矩阵
  s1x=[0 -s1(3) s1(2);s1(3) 0 -s1(1);-s1(2) s1(1) 0];
  s2x=[0 -s2(3) s2(2);s2(3) 0 -s2(1);-s2(2) s2(1) 0];
  s3x=[0 -s3(3) s3(2);s3(3) 0 -s3(1);-s3(2) s3(1) 0];
  s4x=[0 -s4(3) s4(2);s4(3) 0 -s4(1);-s4(2) s4(1) 0];
  s5x=[0 -s5(3) s5(2);s5(3) 0 -s5(1);-s5(2) s5(1) 0];
  s6x=[0 -s6(3) s6(2);s6(3) 0 -s6(1);-s6(2) s6(1) 0];
  % 上铰链位置的反对称矩阵-转动时用
b1x=[0 -br1(3) br1(2);br1(3) 0 -br1(1);-br1(2) br1(1) 0];
b2x=[0 -br2(3) br2(2);br2(3) 0 -br2(1);-br2(2) br2(1) 0];
b3x=[0 -br3(3) br3(2);br3(3) 0 -br3(1);-br3(2) br3(1) 0];
b4x=[0 -br4(3) br4(2);br4(3) 0 -br4(1);-br4(2) br4(1) 0];
b5x=[0 -br5(3) br5(2);br5(3) 0 -br5(1);-br5(2) br5(1) 0];
b6x=[0 -br6(3) br6(2);br6(3) 0 -br6(1);-br6(2) br6(1) 0];

  c1=b1x*s1';
  c2=b2x*s2';
  c3=b3x*s3';
  c4=b4x*s4';
  c5=b5x*s5';
  c6=b6x*s6';
  % 雅可比的转置矩阵
  Jt=[s1' s2' s3' s4' s5' s6';c1 c2 c3 c4 c5 c6];
  % 中间雅可比矩阵3x6
  J1=[eye(3) -b1x];
  J2=[eye(3) -b2x];
  J3=[eye(3) -b3x];
  J4=[eye(3) -b4x];
  J5=[eye(3) -b5x];
  J6=[eye(3) -b6x];
  % 机构的X、V、A
  Vx=[Vp w];
%   Ax=[Ap dw];
  %omiga的反对称阵
  omigax=[0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
  dJ1=[zeros(3) -omigax*b1x+b1x*omigax];
  dJ2=[zeros(3) -omigax*b2x+b2x*omigax];
  dJ3=[zeros(3) -omigax*b3x+b3x*omigax];
  dJ4=[zeros(3) -omigax*b4x+b4x*omigax];
  dJ5=[zeros(3) -omigax*b5x+b5x*omigax];
  dJ6=[zeros(3) -omigax*b6x+b6x*omigax];
% 肢体中间广义坐标到机构的主广义坐标（列向量）
  dx1=J1*Vx';
  dx2=J2*Vx';
  dx3=J3*Vx';
  dx4=J4*Vx';
  dx5=J5*Vx';
  dx6=J6*Vx';
  % 系数
  mce1=1/l1^2*(m11*c11^2+m12*c12^2);
  mce2=1/l2^2*(m11*c11^2+m12*c12^2);
  mce3=1/l3^2*(m11*c11^2+m12*c12^2);
  mce4=1/l4^2*(m11*c11^2+m12*c12^2);
  mce5=1/l5^2*(m11*c11^2+m12*c12^2);
  mce6=1/l6^2*(m11*c11^2+m12*c12^2);
  mco1=1/l1*m12*c12-Ixx/l1^2-mce1;
  mco2=1/l2*m12*c12-Ixx/l2^2-mce2;
  mco3=1/l3*m12*c12-Ixx/l3^2-mce3;
  mco4=1/l4*m12*c12-Ixx/l4^2-mce4;
  mco5=1/l5*m12*c12-Ixx/l5^2-mce5;
  mco6=1/l6*m12*c12-Ixx/l6^2-mce6;
  mge1=1/l1*(m11*c11+m12*(l1-c12));
  mge2=1/l2*(m11*c11+m12*(l2-c12));
  mge3=1/l3*(m11*c11+m12*(l3-c12));
  mge4=1/l4*(m11*c11+m12*(l4-c12));
  mge5=1/l5*(m11*c11+m12*(l5-c12));
  mge6=1/l6*(m11*c11+m12*(l6-c12));
  % 肢体的M、C、G
  M1=m12*(s1'*s1)-1/l1^2*Ixx*s1x*s1x-mce1*s1x*s1x;
  M2=m12*(s2'*s2)-1/l2^2*Ixx*s2x*s2x-mce2*s2x*s2x;
  M3=m12*(s3'*s3)-1/l3^2*Ixx*s3x*s3x-mce3*s3x*s3x;
  M4=m12*(s4'*s4)-1/l4^2*Ixx*s4x*s4x-mce4*s4x*s4x;
  M5=m12*(s5'*s5)-1/l5^2*Ixx*s5x*s5x-mce5*s5x*s5x;
  M6=m12*(s6'*s6)-1/l6^2*Ixx*s6x*s6x-mce6*s6x*s6x;
  C1=-2/l1*mco1*s1*dx1*s1x*s1x+1/l1^2*m12*c12*s1'*(s1x*dx1)'*s1x;
  C2=-2/l2*mco2*s2*dx2*s2x*s2x+1/l2^2*m12*c12*s2'*(s2x*dx2)'*s2x;
  C3=-2/l3*mco3*s3*dx3*s3x*s3x+1/l3^2*m12*c12*s3'*(s3x*dx3)'*s3x;
  C4=-2/l4*mco4*s4*dx4*s4x*s4x+1/l4^2*m12*c12*s4'*(s4x*dx4)'*s4x;
  C5=-2/l5*mco5*s5*dx5*s5x*s5x+1/l5^2*m12*c12*s5'*(s5x*dx5)'*s5x;
  C6=-2/l6*mco6*s6*dx6*s6x*s6x+1/l6^2*m12*c12*s6'*(s6x*dx6)'*s6x;
  G1=(mge1*s1x*s1x-m12*(s1'*s1))*g';
  G2=(mge2*s2x*s2x-m12*(s2'*s2))*g';
  G3=(mge3*s3x*s3x-m12*(s3'*s3))*g';
  G4=(mge4*s4x*s4x-m12*(s4'*s4))*g';
  G5=(mge5*s5x*s5x-m12*(s5'*s5))*g';
  G6=(mge6*s6x*s6x-m12*(s6'*s6))*g';
  % 肢体到机构的M、C、G
Ml1=J1'*M1*J1;
Ml2=J2'*M2*J2;
Ml3=J3'*M3*J3;
Ml4=J4'*M4*J4;
Ml5=J5'*M5*J5;
Ml6=J6'*M6*J6;

Cl1=J1'*M1*dJ1+J1'*C1*J1;
Cl2=J2'*M2*dJ2+J2'*C2*J2;
Cl3=J3'*M3*dJ3+J3'*C3*J3;
Cl4=J4'*M4*dJ4+J4'*C4*J4;
Cl5=J5'*M5*dJ5+J5'*C5*J5;
Cl6=J6'*M6*dJ6+J6'*C6*J6;

Gl1=J1'*G1;
Gl2=J2'*G2;
Gl3=J3'*G3;
Gl4=J4'*G4;
Gl5=J5'*G5;
Gl6=J6'*G6;
% 上平台的M、C、G，
Mp=[m*eye(3) zeros(3);zeros(3) AIp];
Cp=[zeros(3) zeros(3);zeros(3) omigax*AIp];
Gp=[-m*g 0 0 0];
% 机构的M、C、G
M=Mp+Ml1+Ml2+Ml3+Ml4+Ml5+Ml6;
C=Cp+Cl1+Cl2+Cl3+Cl4+Cl5+Cl6;
G=Gp'+Gl1+Gl2+Gl3+Gl4+Gl5+Gl6;
%求解加速度量
F=[y11,y12,y13,y14,y15,y16];
Ax=pinv(M)*(Jt*F'-C*Vx'-G);


end
function Ax1=ax1(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16)
        Ax=ax(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16);
Ax1=Ax(1);
end
function Ax2=ax2(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16)
            Ax=ax(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16);
Ax2=Ax(2);
end
function Ax3=ax3(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16)
        Ax=ax(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16);
Ax3=Ax(3);
end
%计算sx的函数
function Ax4=ax4(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16)
            Ax=ax(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16);
a1=Ax(4);a2=Ax(5);a3=Ax(6);
%旋转矩阵
sa=sin(alpha);
ca=cos(alpha);
sb=sin(beta);
cb=cos(beta);
sg=sin(gama);
cg=cos(gama);
R=[cb*cg sa*sb*cg-ca*sg ca*sb*cg+sg*sa;cb*sg sa*sb*sg+ca*cg ca*sb*sg-cg*sa;-sb sa*cb ca*cb];
%角运动
tr=R(1,1)+R(2,2)+R(3,3);
mid=(tr-1)/2;
theta=acos(mid);
sx=(R(3,2)-R(2,3))/(2*sin(theta));
sy=(R(1,3)-R(3,1))/(2*sin(theta));
sz=(R(2,1)-R(1,2))/(2*sin(theta));
r1d=-bed*sb*cg-gad*sg*cb;
r2d=ald*ca*sb*sg+bed*sa*cb*sg+gad*sa*sb*sg;
r3d=-ald*sa*cg-gad*sg*ca;
r4d=-ald*sa*cb-bed*sb*ca;
trd=r1d+r2d+r3d+r4d;
wd1=0.5*(3-tr^2+2*tr)^(-1.5)*(-2*tr*trd+2*trd)*trd;
wd2=(3-tr^2+2*tr)^(-0.5);
% r1dd=-(x2*sb*cg+bed*bed*cb*cg-bed*sb*gad*sg)-(x3*sg*cb+gad*gad*cg*cb-gad*sg*bed*cb);
% r2dd=x1*ca*sb*sg-ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg
% +x2*sa*cb*sg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg
% +x3*sa*sb*sg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg;
% r3dd=-(x1*sa*cg+ald*ald*ca*cg-ald*sa*gad*sg)-(x3*sg*ca+gad*gad*cg*ca-gad*sg*ald*ca);
% r4dd=-(x1*sa*cb+ald*ald*ca*cb-ald*sa*bed*sb)-(x2*sb*ca+bed*bed*cb*ca-bed*sb*ald*ca);
% trdd=x1*(ca*sb*sg-sa*cg-sa*cb)+x2*(sa*cb*sg-sb*ca-sb*cg)+x3*(sa*sb*sg-sg*cb-sg*ca)-(bed*bed*cb*cg-bed*sb*gad*sg+gad*gad*cg*cb-gad*sg*bed*cb)
% -ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg
% -(ald*ald*ca*cg-ald*sa*gad*sg+gad*gad*cg*ca-gad*sg*ald*ca)-(ald*ald*ca*cb-ald*sa*bed*sb+bed*bed*cb*ca-bed*sb*ald*ca);

m1=ca*sb*sg-sa*cg-sa*cb;
m2=sa*cb*sg-sb*ca-sb*cg;
m3=sa*sb*sg-sg*cb-sg*ca;
m4=-(bed*bed*cb*cg-bed*sb*gad*sg+gad*gad*cg*cb-gad*sg*bed*cb)-ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg-(ald*ald*ca*cg-ald*sa*gad*sg+gad*gad*cg*ca-gad*sg*ald*ca)-(ald*ald*ca*cb-ald*sa*bed*sb+bed*bed*cb*ca-bed*sb*ald*ca);
fun = @(x)root3d(x,sx,sy,sz,wd1,wd2,m1,m2,m3,m4,a1,a2,a3);
x0=[0,0,0];
x=fsolve(fun,x0);
Ax4=x(1);
function F = root3d(x,sx,sy,sz,wd1,wd2,m1,m2,m3,m4,a1,a2,a3)

F(1) = sx*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a1;
F(2) = sy*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a2;
F(3) = sz*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a3;
end
end
%计算sy的函数
function Ax5=ax5(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16)
            Ax=ax(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16);
a1=Ax(4);a2=Ax(5);a3=Ax(6);
%旋转矩阵
sa=sin(alpha);
ca=cos(alpha);
sb=sin(beta);
cb=cos(beta);
sg=sin(gama);
cg=cos(gama);
R=[cb*cg sa*sb*cg-ca*sg ca*sb*cg+sg*sa;cb*sg sa*sb*sg+ca*cg ca*sb*sg-cg*sa;-sb sa*cb ca*cb];
%角运动
tr=R(1,1)+R(2,2)+R(3,3);
mid=(tr-1)/2;
theta=acos(mid);
sx=(R(3,2)-R(2,3))/(2*sin(theta));
sy=(R(1,3)-R(3,1))/(2*sin(theta));
sz=(R(2,1)-R(1,2))/(2*sin(theta));
r1d=-bed*sb*cg-gad*sg*cb;
r2d=ald*ca*sb*sg+bed*sa*cb*sg+gad*sa*sb*sg;
r3d=-ald*sa*cg-gad*sg*ca;
r4d=-ald*sa*cb-bed*sb*ca;
trd=r1d+r2d+r3d+r4d;
wd1=0.5*(3-tr^2+2*tr)^(-1.5)*(-2*tr*trd+2*trd)*trd;
wd2=(3-tr^2+2*tr)^(-0.5);
% r1dd=-(x2*sb*cg+bed*bed*cb*cg-bed*sb*gad*sg)-(x3*sg*cb+gad*gad*cg*cb-gad*sg*bed*cb);
% r2dd=x1*ca*sb*sg-ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg
% +x2*sa*cb*sg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg
% +x3*sa*sb*sg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg;
% r3dd=-(x1*sa*cg+ald*ald*ca*cg-ald*sa*gad*sg)-(x3*sg*ca+gad*gad*cg*ca-gad*sg*ald*ca);
% r4dd=-(x1*sa*cb+ald*ald*ca*cb-ald*sa*bed*sb)-(x2*sb*ca+bed*bed*cb*ca-bed*sb*ald*ca);
% trdd=x1*(ca*sb*sg-sa*cg-sa*cb)+x2*(sa*cb*sg-sb*ca-sb*cg)+x3*(sa*sb*sg-sg*cb-sg*ca)-(bed*bed*cb*cg-bed*sb*gad*sg+gad*gad*cg*cb-gad*sg*bed*cb)
% -ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg
% -(ald*ald*ca*cg-ald*sa*gad*sg+gad*gad*cg*ca-gad*sg*ald*ca)-(ald*ald*ca*cb-ald*sa*bed*sb+bed*bed*cb*ca-bed*sb*ald*ca);

m1=ca*sb*sg-sa*cg-sa*cb;
m2=sa*cb*sg-sb*ca-sb*cg;
m3=sa*sb*sg-sg*cb-sg*ca;
m4=-(bed*bed*cb*cg-bed*sb*gad*sg+gad*gad*cg*cb-gad*sg*bed*cb)-ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg-(ald*ald*ca*cg-ald*sa*gad*sg+gad*gad*cg*ca-gad*sg*ald*ca)-(ald*ald*ca*cb-ald*sa*bed*sb+bed*bed*cb*ca-bed*sb*ald*ca);
fun = @(x)root3d(x,sx,sy,sz,wd1,wd2,m1,m2,m3,m4,a1,a2,a3);
x0=[0,0,0];
x=fsolve(fun,x0);
Ax5=x(2);
function F = root3d(x,sx,sy,sz,wd1,wd2,m1,m2,m3,m4,a1,a2,a3)

F(1) = sx*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a1;
F(2) = sy*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a2;
F(3) = sz*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a3;
end
end
%计算sz的函数
function Ax6=ax6(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16)
            Ax=ax(y1,y2,y3,alpha,beta,gama,dy1,dy2,dy3,ald,bed,gad,y11,y12,y13,y14,y15,y16);
a1=Ax(4);a2=Ax(5);a3=Ax(6);
%旋转矩阵
sa=sin(alpha);
ca=cos(alpha);
sb=sin(beta);
cb=cos(beta);
sg=sin(gama);
cg=cos(gama);
R=[cb*cg sa*sb*cg-ca*sg ca*sb*cg+sg*sa;cb*sg sa*sb*sg+ca*cg ca*sb*sg-cg*sa;-sb sa*cb ca*cb];
%角运动
tr=R(1,1)+R(2,2)+R(3,3);
mid=(tr-1)/2;
theta=acos(mid);
sx=(R(3,2)-R(2,3))/(2*sin(theta));
sy=(R(1,3)-R(3,1))/(2*sin(theta));
sz=(R(2,1)-R(1,2))/(2*sin(theta));
r1d=-bed*sb*cg-gad*sg*cb;
r2d=ald*ca*sb*sg+bed*sa*cb*sg+gad*sa*sb*sg;
r3d=-ald*sa*cg-gad*sg*ca;
r4d=-ald*sa*cb-bed*sb*ca;
trd=r1d+r2d+r3d+r4d;
wd1=0.5*(3-tr^2+2*tr)^(-1.5)*(-2*tr*trd+2*trd)*trd;
wd2=(3-tr^2+2*tr)^(-0.5);
% r1dd=-(x2*sb*cg+bed*bed*cb*cg-bed*sb*gad*sg)-(x3*sg*cb+gad*gad*cg*cb-gad*sg*bed*cb);
% r2dd=x1*ca*sb*sg-ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg
% +x2*sa*cb*sg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg
% +x3*sa*sb*sg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg;
% r3dd=-(x1*sa*cg+ald*ald*ca*cg-ald*sa*gad*sg)-(x3*sg*ca+gad*gad*cg*ca-gad*sg*ald*ca);
% r4dd=-(x1*sa*cb+ald*ald*ca*cb-ald*sa*bed*sb)-(x2*sb*ca+bed*bed*cb*ca-bed*sb*ald*ca);
% trdd=x1*(ca*sb*sg-sa*cg-sa*cb)+x2*(sa*cb*sg-sb*ca-sb*cg)+x3*(sa*sb*sg-sg*cb-sg*ca)-(bed*bed*cb*cg-bed*sb*gad*sg+gad*gad*cg*cb-gad*sg*bed*cb)
% -ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg
% -(ald*ald*ca*cg-ald*sa*gad*sg+gad*gad*cg*ca-gad*sg*ald*ca)-(ald*ald*ca*cb-ald*sa*bed*sb+bed*bed*cb*ca-bed*sb*ald*ca);

m1=ca*sb*sg-sa*cg-sa*cb;
m2=sa*cb*sg-sb*ca-sb*cg;
m3=sa*sb*sg-sg*cb-sg*ca;
m4=-(bed*bed*cb*cg-bed*sb*gad*sg+gad*gad*cg*cb-gad*sg*bed*cb)-ald*ald*sa*sb*sg+bed*ald*ca*cb*sg+gad*ald*ca*sb*cg+bed*ald*ca*cb*sg-bed*sa*bed*sb*sg+bed*sa*cb*gad*cg+gad*ald*ca*sb*sg+gad*sa*bed*cb*sg+gad*sa*sb*gad*cg-(ald*ald*ca*cg-ald*sa*gad*sg+gad*gad*cg*ca-gad*sg*ald*ca)-(ald*ald*ca*cb-ald*sa*bed*sb+bed*bed*cb*ca-bed*sb*ald*ca);
fun = @(x)root3d(x,sx,sy,sz,wd1,wd2,m1,m2,m3,m4,a1,a2,a3);
x0=[0,0,0];
x=fsolve(fun,x0);
Ax6=x(3);
function F = root3d(x,sx,sy,sz,wd1,wd2,m1,m2,m3,m4,a1,a2,a3)

F(1) = sx*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a1;
F(2) = sy*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a2;
F(3) = sz*(wd1-wd2*(m1*x(1)+m2*x(2)+m3*x(3)+m4))-a3;
end
end
