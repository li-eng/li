
% 上平台铰链位置
tic

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
m_dj=0.667417;m_sg=0.072068;m_gz=0.07012;
m11=m_dj+2*m_gz+m_sg;m12=2*m_gz;
c11=0.029811+0.024+0.0365;c12=0.09+0.021;%+0.046
g=[0 0 -9.8];
% 上平台质量及惯性矩
m=0.8151490458;Ip=[0.00403852037 0 0;0 0.0040385 0;0 0 0.0080758];
% 上平台中心位置
t=0:0.0001:1;
% z=0.3825+0.09*t.^2-0.06*t.^3;%z向移动30mm
% z1=0.18*t-0.18*t.^2;
% z2=0.18-0.36*t;
%y向移动
% y=0.09*t.^2-0.06*t.^3;%y向移动30mm
% y1=0.18*t-0.18*t.^2;
% y2=0.18-0.36*t;
%x向移动30mm
% x=0.09*t.^2-0.06*t.^3;
% x1=0.18*t-0.18*t.^2;
% x2=0.18-0.36*t;
% theta=pi/2*t.^2-pi/3*t.^3;%绕z-1s旋转30度
% omiga=pi*t-pi*t.^2;
% omigado=pi-2*pi*t;
theta=pi/6*t.^2-pi/9*t.^3;%绕z-1s旋转10度
omiga=pi/3*t-pi/3*t.^2;
omigado=pi/3-2*pi/3*t;
sx=1;sy=1;sz=0;
N=length(t);
F=zeros(6,N);

% theta=zeros(1,N);omiga=zeros(1,N);omigado=zeros(1,N);

for i=1:N
%     p=[x(i) 0 0.3825];Vp=[x1(i) 0 0];Ap=[x2(i) 0 0];%沿x直线运动
%    p=[0 y(i) 0.3825];Vp=[0 y1(i) 0];Ap=[0 y2(i) 0];%沿y直线运动
%  p=[0 0 z(i)];Vp=[0 0 z1(i)];Ap=[0 0 z2(i)];%沿z直线运动
p=[0 0 0.3825];Vp=[0 0 0];Ap=[0 0 0];%旋转时用

% theta(i)=0;omiga(i)=0;omigado(i)=0;%直线运动用
sth=sin(theta(i));
cth=cos(theta(i));
vth=1-cth;
%旋转矩阵
R=[sx*sx*vth+cth sx*sy*vth-sz*sth sx*sz*vth+sy*sth;sy*sx*vth+sz*sth sy*sy*vth+cth sy*sz*vth-sx*sth;sz*sx*vth-sy*sth sz*sy*vth+sx*sth sz*sz*vth+cth];

%角运动
omiga1=omiga(i)*sx;
omiga2=omiga(i)*sy;
omiga3=omiga(i)*sz;
omigad1=omigado(i)*sx;
omigad2=omigado(i)*sy;
omigad3=omigado(i)*sz;
 %静坐标系下的上铰链点位置
 br1=R*b1'; br2=R*b2'; br3=R*b3'; br4=R*b4'; br5=R*b5'; br6=R*b6';%列向量
   %杆的转动惯量
  Ixx=0.003395+0.000391;
  % 杆长及其单位向量
  l1=sqrt(p*p'+br1'*br1+a1*a1'-2*p*a1'+2*p*br1-2*br1'*a1');
  l2=sqrt(p*p'+br2'*br2+a2*a2'-2*p*a2'+2*p*br2-2*br2'*a2');
  l3=sqrt(p*p'+br3'*br3+a3*a3'-2*p*a3'+2*p*br3-2*br3'*a3');
  l4=sqrt(p*p'+br4'*br4+a4*a4'-2*p*a4'+2*p*br4-2*br4'*a4');
  l5=sqrt(p*p'+br5'*br5+a5*a5'-2*p*a5'+2*p*br5-2*br5'*a5');
  l6=sqrt(p*p'+br6'*br6+a6*a6'-2*p*a6'+2*p*br6-2*br6'*a6');

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

% 上铰链位置的反对称矩阵-无旋转时用
%   b1x=[0 0.0365 -0.07320508076;-0.0365 0 -0.07320508076;0.07320508076 0.07320508076 0];
%   b2x=[0 0.0365 -0.07320508076;-0.0365 0 0.07320508076;0.07320508076 -0.07320508076 0];
%   b3x=[0 0.0365 -0.02679491924;-0.0365 0 0.1;0.02679491924 -0.100 0];
%   b4x=[0 0.0365 0.100;-0.0365 0 0.02679491924;-0.100 -0.02679491924 0];
%   b5x=[0 0.0365 0.100;-0.0365 0 -0.02679491924;-0.100 0.02679491924 0];
%   b6x=[0 0.0365 -0.02679491924;-0.0365 0 -0.100;0.02679491924 0.100 0];
% 上铰链位置的反对称矩阵
b1x=[0 -br1(3) br1(2);br1(3) 0 -br1(1);-br1(2) br1(1) 0];
b2x=[0 -br2(3) br2(2);br2(3) 0 -br2(1);-br2(2) br2(1) 0];
b3x=[0 -br3(3) br3(2);br3(3) 0 -br3(1);-br3(2) br3(1) 0];
b4x=[0 -br4(3) br4(2);br4(3) 0 -br4(1);-br4(2) br4(1) 0];
b5x=[0 -br5(3) br5(2);br5(3) 0 -br5(1);-br5(2) br5(1) 0];
b6x=[0 -br6(3) br6(2);br6(3) 0 -br6(1);-br6(2) br6(1) 0];
% 矩阵叉乘
  c1=b1x*s1';
  c2=b2x*s2';
  c3=b3x*s3';
  c4=b4x*s4';
  c5=b5x*s5';
  c6=b6x*s6';

% 机构的X、V、A
w=[omiga1 omiga2 omiga3];
 Vx=[Vp w];Ax=[Ap omigad1 omigad2 omigad3];
% 雅可比的转置矩阵
  Jt=[s1' s2' s3' s4' s5' s6';c1 c2 c3 c4 c5 c6];

% 中间雅可比矩阵3x6
  J1=[eye(3) -b1x];
  J2=[eye(3) -b2x];
  J3=[eye(3) -b3x];
  J4=[eye(3) -b4x];
  J5=[eye(3) -b5x];
  J6=[eye(3) -b6x];
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

% 上平台的M、C、G，其中C为0
Mp=[m*eye(3) zeros(3);zeros(3) Ip];
% Cp=[zeros(3) zeros(3);zeros(3) ];
Gp=[-m*g 0 0 0];
% 机构的M、C、G
M=Mp+Ml1+Ml2+Ml3+Ml4+Ml5+Ml6;
C=Cl1+Cl2+Cl3+Cl4+Cl5+Cl6;
G=Gp'+Gl1+Gl2+Gl3+Gl4+Gl5+Gl6;

f=Jt\(M*Ax'+C*Vx'+G);

F(:,i)=f;  
end
xlswrite('f_rxy10d_0.0001.xlsx',F);
 plot(t,F(1,:),'r-');
 hold on
 plot(t,F(2,:),'b-.');
 hold on 
 plot(t,F(3,:),'g--');
 hold on
 plot(t,F(4,:),'k*');
 hold on
 plot(t,F(5,:),'y');
  hold on
 plot(t,F(6,:),'+');
  legend('l1-','l2-.','l3--','l4*','l5y','l6+')

%  A1 = load('force1z_1s_30mm.tab');
%  A2 = load('force_5mm.tab');
%  A3 = load('force_5mm.tab');
%  A4 = load('force_5mm.tab');
%  A5 = load('force_5mm.tab');
%  A6 = load('force_5mm.tab');
% B=load('force_rotationz30d_m01.tab');
% C=load('rotationz30d_02.tab');
% C1=load('rotationz30d_03.tab');
% C2=load('y5mm_04.tab');
%   plot(t,A1(:,2));
% % hold on
% error01=A1(:,2)-F(6,:)';
% error02=A2(:,2)-F(2,:)';
% error03=A3(:,2)-F(3,:)';
% error04=A4(:,2)-F(4,:)';
% error05=A5(:,2)-F(5,:)';
% error06=A6(:,2)-F(6,:)';
%   error=B(45:end,2)-F(1,45:end)';
%   error1=C1(:,2)-F(3,:)';
%   error2=C2(:,2)-F(4,:)';
%  u=0:0.01:0.56;
%   plot(t,error01)
%   plot(t,error01,'--',t,error02,'-.',t,error03,'*',t,error04,t,error05,t,error06,'o');
%   xlabel('t/s');      %写出x坐标标示
%   ylabel('force-error');      %写出y坐标标示
%   title('rotation-z-30d-m01');  %写出此图形代表标题
%   legend('rod-1','rod-2','rod-3','rod-4','rod-5','rod-6')
% legend('rod-1、3、5','rod-2、4、6')
toc