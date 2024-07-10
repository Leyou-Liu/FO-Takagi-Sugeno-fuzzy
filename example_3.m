clear;
clc;
%%%              
A11=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 17231/16951 -2/9 0;0 0 0 0 0 0];
A12=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 1577/16951 -2/9 0;0 0 0 0 0 0];
Ad=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 -7/45;0 0 0 0 0 -10];
% B=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
% B=[0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];%%%             
%% 
B=eye(6);  %% 
I=eye(6);
setlmis([]);
Q=lmivar(1,[6 1]);
M1=lmivar(2,[6 6]);
M2=lmivar(2,[6 6]);
S=lmivar(1,[6 1]);
lmiterm([1 1 1 Q],1,A11','s');
lmiterm([1 1 1 -M1],1,B','s');
lmiterm([1 1 2 0],Ad);
lmiterm([1 1 3 Q],1,1);

lmiterm([1 2 2 S],-1,1);
lmiterm([1 2 3 0],0);

lmiterm([1 3 3 S],1,1);
lmiterm([1 3 3 0],-2*I);

lmiterm([2 1 1 Q],1,A12','s');
lmiterm([2 1 1 -M2],1,B','s');
lmiterm([2 1 2 0],Ad);
lmiterm([2 1 3 Q],1,1);

lmiterm([2 2 2 S],-1,1);
lmiterm([2 2 3 0],0);

lmiterm([2 3 3 S],1,1);
lmiterm([2 3 3 0],-2*I);

lmiterm([3 1 1 Q],1,A11'+A12','s');
lmiterm([3 1 1 -M1],1,B','s');
lmiterm([3 1 1 -M2],1,B','s');

lmiterm([3 1 2 0],2*Ad);
lmiterm([3 1 3 Q],2,1);

lmiterm([3 2 2 S],-2,1);
lmiterm([3 2 3 0],0);

lmiterm([3 3 3 S],2,1);
lmiterm([3 3 3 0],-4*I);

lmiterm([-4 1 1 Q],1,1);
lmiterm([-5 1 1 S],1,1);
lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
MATQ=dec2mat(lmisys,xfeas,Q);
MATS=dec2mat(lmisys,xfeas,S);
MATM1=dec2mat(lmisys,xfeas,M1);
MATM2=dec2mat(lmisys,xfeas,M2);

muk1=MATM1/MATQ;
muk2=MATM2/MATQ;

Tya=0.1;
A11=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 17231/16951 -2/9 0;0 0 0 0 0 0];
A12=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 1577/16951 -2/9 0;0 0 0 0 0 0];
Ad=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 -7/45;0 0 0 0 0 -1/Tya];
B=eye(6);
% B=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];%g
% B=[0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
I=eye(6);C=eye(6);
% Bw=1.04*eye(6);
Bw=1*eye(6);yip=0.08;ga=1;   %0.2 xt
setlmis([]);
Q=lmivar(1,[6 1]);
M1=lmivar(2,[6 6]);
M2=lmivar(2,[6 6]);
S=lmivar(1,[6 1]);
lmiterm([1 1 1 Q],1,A11','s');
lmiterm([1 1 1 -M1],1,B','s');
lmiterm([1 1 1 S],1,1);
lmiterm([1 1 1 0],C'*C);
lmiterm([1 1 1 0],-I);
lmiterm([1 1 2 Q],1,Ad);
lmiterm([1 1 3 Q],1,Bw);
lmiterm([1 1 4 0],I);
lmiterm([1 2 2 S],-1,1);
lmiterm([1 2 3 0],0);
lmiterm([1 2 4 0],0);
lmiterm([1 3 3 0],-ga*ga*I);
lmiterm([1 3 4 0],0);
lmiterm([1 4 4 0],-I);
lmiterm([2 1 1 Q],1,A12','s');
lmiterm([2 1 1 -M2],1,B','s');
lmiterm([2 1 1 S],1,1);
lmiterm([2 1 1 0],C'*C);
lmiterm([2 1 1 0],-I);
lmiterm([2 1 2 Q],1,Ad);
lmiterm([2 1 3 Q],1,Bw);
lmiterm([2 1 4 0],I);
lmiterm([2 2 2 S],-1,1);
lmiterm([2 2 3 0],0);
lmiterm([2 2 4 0],0);
lmiterm([2 3 3 0],-ga*ga*I);
lmiterm([2 3 4 0],0);
lmiterm([2 4 4 0],-I);
lmiterm([3 1 1 Q],1,A11'+A12','s');
lmiterm([3 1 1 -M1],1,B','s');
lmiterm([3 1 1 -M2],1,B','s');
lmiterm([3 1 1 S],1,1);
lmiterm([3 1 1 0],C'*C);
lmiterm([3 1 1 0],-I);
lmiterm([3 1 2 Q],1,Ad);
lmiterm([3 1 3 Q],1,Bw);
lmiterm([3 1 4 0],I);
lmiterm([3 2 2 S],-1,1);
lmiterm([3 2 3 0],0);
lmiterm([3 2 4 0],0);
lmiterm([3 3 3 0],-ga*ga*I);
lmiterm([3 3 4 0],0);
lmiterm([3 4 4 0],-I);

lmiterm([-4 1 1 Q],1,1);
lmiterm([-5 1 1 S],1,1);
lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
MATQ=dec2mat(lmisys,xfeas,Q);
MATS=dec2mat(lmisys,xfeas,S);
MATM1=dec2mat(lmisys,xfeas,M1);
MATM2=dec2mat(lmisys,xfeas,M2);
tmin;
k1=MATM1/MATQ;
k2=MATM2/MATQ;

LL1 =[-9.5283  -29.3691    9.3871   18.5425   75.2911   -3.0387
   28.3691   -9.5283   10.5806    0.9130  -35.0090 -185.2046
   14.8506   28.2480   -3.3178    2.9739 -107.2577  -36.8286];
LL2 =[-9.5283   39.6651    8.9416   -1.6484   12.8088  -20.3104
  -40.6651   -9.5283   36.0380   14.1962    7.0278 -152.7111
   16.8373    0.1991   -8.0356   -0.3294  -18.2063   17.0228];
kk1 =[0.0094    0.0034   -0.2119   -9.3270 -157.5028    0.0603
   -6.4301   -0.0310   -0.7939 -157.4898   -9.1058    0.0349
   -2.0080   -2.0222    1.9866    1.4774   -0.3522  -25.8593];
kk2 =[0.0095    0.0034   -0.2136   -9.3270 -157.0453    0.0732
   -6.4297   -0.0309   -0.8043 -157.0236   -9.1047    0.1578
   -2.0080   -2.0222    1.9867    1.4769   -0.3993  -25.8605];
BB=[0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
zhenggti=BB*(kk1+kk2);

% k1=[0 0 0 0 0 0
%    0 0 0 0 0 0
%    0 0 0 0 0 0
%   0.0094    0.0034   -0.2119   -9.3270 -157.5028    0.0603
%    -6.4301   -0.0310   -0.7939 -157.4898   -9.1058    0.0349
%    -2.0080   -2.0222    1.9866    1.4774   -0.3522  -25.8593];
% 
% k2=[0 0 0 0 0 0
%    0 0 0 0 0 0
%    0 0 0 0 0 0
%     0.0095    0.0034   -0.2136   -9.3270 -157.0453    0.0732
%    -6.4297   -0.0309   -0.8043 -157.0236   -9.1047    0.1578
%    -2.0080   -2.0222    1.9867    1.4769   -0.3993  -25.8605]; %

KKK1 =[0.0094    0.0034   -0.2119   -9.3270 -157.5028    0.0603
   -6.4301   -0.0310   -0.7939 -157.4898   -9.1058    0.0349
   -2.0080   -2.0222    1.9866    1.4774   -0.3522  -25.8593];
KKK2 =[0.0095    0.0034   -0.2136   -9.3270 -157.0453    0.0732
   -6.4297   -0.0309   -0.8043 -157.0236   -9.1047    0.1578
   -2.0080   -2.0222    1.9867    1.4769   -0.3993  -25.8605];

zhengtii1=diag(0.6*(muk1)+0.4*(muk2));
zhengti1=diag(zhengtii1);
v1=zhengti1(1,1); v2=zhengti1(1,2); v3=zhengti1(1,3); v4=zhengti1(1,4);v5=zhengti1(1,5);v6=zhengti1(1,6);
v7=zhengti1(2,1); v8=zhengti1(2,2); v9=zhengti1(2,3); v10=zhengti1(2,4);v11=zhengti1(2,5);v12=zhengti1(2,6);
v13=zhengti1(3,1);v14=zhengti1(3,2);v15=zhengti1(3,3);v16=zhengti1(3,4);v17=zhengti1(3,5);v18=zhengti1(3,6);
v19=zhengti1(4,1);v20=zhengti1(4,2);v21=zhengti1(4,3);v22=zhengti1(4,4);v23=zhengti1(4,5);v24=zhengti1(4,6);
v25=zhengti1(5,1);v26=zhengti1(5,2);v27=zhengti1(5,3);v28=zhengti1(5,4);v29=zhengti1(5,5);v30=zhengti1(5,6);
v31=zhengti1(6,1);v32=zhengti1(6,2);v33=zhengti1(6,3);v34=zhengti1(6,4);v35=zhengti1(6,5);v36=zhengti1(6,6);

q1=0.98;q2=0.98;q3=0.98;q4=0.98;q5=0.98;q6=0.98;  
% q1=0.9;q2=0.9;q3=0.9;q4=0.9;q5=0.9;q6=0.9;  

Tya=0.1; td1=0.1;

h=0.01;N=1000;           
x0=0.01;y0=0.01;z0=0.01;w0=0.01;g0=0.01;m0=0.01;   %

x(N+1)=[0];y(N+1)=[0];z(N+1)=[0];w(N+1)=[0];g(N+1)=[0];m(N+1)=[0];  %
x1(N+1)=[0];y1(N+1)=[0];z1(N+1)=[0];w1(N+1)=[0];g1(N+1)=[0];m1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 )
x1(1)=x0+h^q1*(  y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0   )/(gamma(q1)*q1);
y1(1)=y0+h^q2*(  z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0   )/(gamma(q2)*q2);
z1(1)=z0+h^q3*(  -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0    )/(gamma(q3)*q3);
w1(1)=w0+h^q4*(  314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   )/(gamma(q4)*q4);
g1(1)=g0+h^q5*(  6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0    +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0   )/(gamma(q5)*q5);
m1(1)=m0+h^q6*(  -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0   )/(gamma(q6)*q6);
%    f( x1(1) , y1(1) , z1(1)      f( x0 , y0 , z0 
x(1)=x0+h^q1*((  y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0  )+q1*(  y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0    ))/gamma(q1+2);
y(1)=y0+h^q2*((    z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0   )+q2*(  z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0     ))/gamma(q2+2);
z(1)=z0+h^q3*((    -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0   )+q3*(  -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0   ))/gamma(q3+2);
w(1)=w0+h^q4*((   314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   )+q4*(  314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   ))/gamma(q4+2);
g(1)=g0+h^q5*((   6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0     )+q5*(  6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0   ))/gamma(q5+2);
m(1)=m0+h^q6*((   -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0  )+q6*(  -10*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0   ))/gamma(q6+2);
%++++++++++++++++++++++++ ++++++
for n=0:td1/h-1               %    f( x0 , y0 , z0 
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(     y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0     );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0  );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(    6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0   );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(   -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0    );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(    y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0  );
N2=((n+1)^q2-n^q2)*(     z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0     );
N3=((n+1)^q3-n^q3)*(    -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0 );
N4=((n+1)^q4-n^q4)*(      314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0     );
N5=((n+1)^q5-n^q5)*(   6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0   );
N6=((n+1)^q6-n^q6)*(     -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0       );
for j=1:n   %    f(  x(j) , y(j) , z(j) 
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)   );
M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)     );
M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*x(j)-24*y(j)-3*z(j)+m(j)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)   );
M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)   );
M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(    6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)      );
M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(  -1/Tya*m(j)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)    );
      %    f(  x(j) , y(j) , z(j)  
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(      y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)    );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(    z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)  );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(    -24*x(j)-24*y(j)-3*z(j)+m(j)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)    );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(     314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)    );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(    6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)     );
        N6=N6+((n-j+1)^q6-(n-j)^q6)*(    -1/Tya*m(j)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)      );
end   
x1(n+1)=x0+h^q1*N1/(gamma(q1)*q1);                        %2）
y1(n+1)=y0+h^q2*N2/(gamma(q2)*q2);
z1(n+1)=z0+h^q3*N3/(gamma(q3)*q3);
w1(n+1)=w0+h^q4*N4/(gamma(q4)*q4);
g1(n+1)=g0+h^q5*N5/(gamma(q5)*q5);
m1(n+1)=m0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  
x(n+2)=x0+h^q1*(   y(n+1)      +v1*x(n+1)+v2*y(n+1)+v3*z(n+1)+v4*w(n+1)+v5*g(n+1)+v6*m(n+1)   +M1)/gamma(q1+2);
y(n+2)=y0+h^q2*(    z(n+1)      +v7*x(n+1)+v8*y(n+1)+v9*z(n+1)+v10*w(n+1)+v11*g(n+1)+v12*m(n+1)        +M2)/gamma(q2+2);
z(n+2)=z0+h^q3*(     -24*x(n+1)-24*y(n+1)-3*z(n+1)+m(n+1)        +v13*x(n+1)+v14*y(n+1)+v15*z(n+1)+v16*w(n+1)+v17*g(n+1)+v18*m(n+1)    +M3)/gamma(q3+2);
w(n+2)=w0+h^q4*(     314*g(n+1)    +v19*x(n+1)+v20*y(n+1)+v21*z(n+1)+v22*w(n+1)+v23*g(n+1)+v24*m(n+1)        +M4)/gamma(q4+2);
g(n+2)=g0+h^q5*(     6.4*x(n+1)+0.8*z(n+1)-3/23*sin(w(n+1))+320/30429*sin(2*w(n+1))-2/9*g(n+1)-7/45*m(n+1)  +v25*x(n+1)+v26*y(n+1)+v27*z(n+1)+v28*w(n+1)+v29*g(n+1)+v30*m(n+1)    +M5)/gamma(q5+2);
m(n+2)=m0+h^q6*(     -1/Tya*m(n+1)   +v31*x(n+1)+v32*y(n+1)+v33*z(n+1)+v34*w(n+1)+v35*g(n+1)+v36*m(n+1)      +M6)/gamma(q6+2);

end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%------------------
for n=td1/h:N
    %                    f( x0 , y0 , z0 
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(    y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(      z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0     );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0   );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(     6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0  );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(       -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0    );

N1=((n+1)^q1-n^q1)*(     y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0       );
N2=((n+1)^q2-n^q2)*(      z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0      );
N3=((n+1)^q3-n^q3)*(     -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0  );
N4=((n+1)^q4-n^q4)*(    314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   );
N5=((n+1)^q5-n^q5)*(     6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0  );
N6=((n+1)^q6-n^q6)*(     -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0  );
for j=1:td1/h      
       %                                       f(  x(j) , y(j) , z(j)   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(      y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)      );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)  );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -24*x(j)-24*y(j)-3*z(j)+m(j)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)  );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(       314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)     );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(     6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)  );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(       -1/Tya*m(j)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)     );
      %    f(  x(j) , y(j) , z(j) 
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*x(j)-24*y(j)-3*z(j)+m(j)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)  );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(       314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)       );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(       6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j) );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(        -1/Tya*m(j)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)       );
end   
for j=td1/h+1:n               
       %                    f(  x(j) , y(j) , z(j)       
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)     );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*x(j)-24*y(j)-3*z(j)+m(j-td1/h)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j) );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(      314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)      );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(       6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j-td1/h)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)   );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(       -1/Tya*m(j-td1/h)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)    );      
      %    f(  x(j) , y(j) , z(j)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)      );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)   );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*x(j)-24*y(j)-3*z(j)+m(j-td1/h)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)   );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(      314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)   );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(      6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j-td1/h)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)  );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(       -1/Tya*m(j-td1/h)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)    );      
end
x1(n+1)=x0+h^q1*N1/(gamma(q1)*q1);                        %2）
y1(n+1)=y0+h^q2*N2/(gamma(q2)*q2);
z1(n+1)=z0+h^q3*N3/(gamma(q2)*q3);
w1(n+1)=w0+h^q4*N4/(gamma(q2)*q4);
g1(n+1)=g0+h^q5*N5/(gamma(q5)*q5);
m1(n+1)=m0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) )
x(n+1)=x0+h^q1*(     y(n+1)      +v1*x(n+1)+v2*y(n+1)+v3*z(n+1)+v4*w(n+1)+v5*g(n+1)+v6*m(n+1)   +M1)/gamma(q1+2);
y(n+1)=y0+h^q2*(     z(n+1)      +v7*x(n+1)+v8*y(n+1)+v9*z(n+1)+v10*w(n+1)+v11*g(n+1)+v12*m(n+1)     +M2)/gamma(q2+2);
z(n+1)=z0+h^q3*(     -24*x(n+1)-24*y(n+1)-3*z(n+1)+m(n+1-td1/h)        +v13*x(n+1)+v14*y(n+1)+v15*z(n+1)+v16*w(n+1)+v17*g(n+1)+v18*m(n+1)    +M3)/gamma(q3+2);
w(n+1)=w0+h^q4*(     314*g(n+1)    +v19*x(n+1)+v20*y(n+1)+v21*z(n+1)+v22*w(n+1)+v23*g(n+1)+v24*m(n+1)      +M4)/gamma(q4+2);
g(n+1)=g0+h^q5*(     6.4*x(n+1)+0.8*z(n+1)-3/23*sin(w(n+1))+320/30429*sin(2*w(n+1))-2/9*g(n+1)-7/45*m(n+1-td1/h)  +v25*x(n+1)+v26*y(n+1)+v27*z(n+1)+v28*w(n+1)+v29*g(n+1)+v30*m(n+1)    +M5)/gamma(q5+2);
m(n+1)=m0+h^q6*(     -1/Tya*m(n+1-td1/h)   +v31*x(n+1)+v32*y(n+1)+v33*z(n+1)+v34*w(n+1)+v35*g(n+1)+v36*m(n+1)      +M6)/gamma(q6+2);
end


zhengtii2=diag(0.6*(muk1)+0.4*(muk2));
zhengti2=diag(zhengtii2);

vv1=zhengti2(1,1); vv2=zhengti2(1,2); vv3=zhengti2(1,3); vv4=zhengti2(1,4);vv5=zhengti2(1,5);vv6=zhengti2(1,6);
vv7=zhengti2(2,1); vv8=zhengti2(2,2); vv9=zhengti2(2,3); vv10=zhengti2(2,4);vv11=zhengti2(2,5);vv12=zhengti2(2,6);
vv13=zhengti2(3,1);vv14=zhengti2(3,2);vv15=zhengti2(3,3);vv16=zhengti2(3,4);vv17=zhengti2(3,5);vv18=zhengti2(3,6);
vv19=zhengti2(4,1);vv20=zhengti2(4,2);vv21=zhengti2(4,3);vv22=zhengti2(4,4);vv23=zhengti2(4,5);vv24=zhengti2(4,6);
vv25=zhengti2(5,1);vv26=zhengti2(5,2);vv27=zhengti2(5,3);vv28=zhengti2(5,4);vv29=zhengti2(5,5);vv30=zhengti2(5,6);
vv31=zhengti2(6,1);vv32=zhengti2(6,2);vv33=zhengti2(6,3);vv34=zhengti2(6,4);vv35=zhengti2(6,5);vv36=zhengti2(6,6);

q1=0.95;q2=0.95;q3=0.95;q4=0.95;q5=0.95;q6=0.95; 
% q1=0.85;q2=0.85;q3=0.85;q4=0.85;q5=0.85;q6=0.85;  

Tya=0.1; td2=0.1;

xx0=0.01;yy0=0.01;zz0=0.01;ww0=0.01;gg0=0.01;mm0=0.01;   %initial value

xx(N+1)=[0];yy(N+1)=[0];zz(N+1)=[0];ww(N+1)=[0];gg(N+1)=[0];mm(N+1)=[0];  %efficiency need improve
xx1(N+1)=[0];yy1(N+1)=[0];zz1(N+1)=[0];ww1(N+1)=[0];gg1(N+1)=[0];mm1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
xx1(1)=xx0+h^q1*(  yy0      +vv1*xx0+vv2*yy0+vv3*zz0+vv4*ww0+vv5*gg0+vv6*mm0   )/(gamma(q1)*q1);
yy1(1)=yy0+h^q2*(  zz0      +vv7*xx0+vv8*yy0+vv9*zz0+vv10*ww0+vv11*gg0+vv12*mm0   )/(gamma(q2)*q2);
zz1(1)=zz0+h^q3*(  -24*xx0-24*yy0-3*zz0+mm0        +vv13*xx0+vv14*yy0+vv15*zz0+vv16*ww0+vv17*gg0+vv18*mm0    )/(gamma(q3)*q3);
ww1(1)=ww0+h^q4*(  314*gg0    +vv19*xx0+vv20*yy0+vv21*zz0+vv22*ww0+vv23*gg0+vv24*mm0   )/(gamma(q4)*q4);
gg1(1)=gg0+h^q5*(  6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0    +vv25*xx0+vv26*yy0+vv27*zz0+vv28*ww0+vv29*gg0+vv30*mm0   )/(gamma(q5)*q5);
mm1(1)=mm0+h^q6*(  -1/Tya*mm0   +vv31*xx0+vv32*yy0+vv33*zz0+vv34*ww0+vv35*gg0+vv36*mm0   )/(gamma(q6)*q6);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
xx(1)=xx0+h^q1*((  yy0      +vv1*xx0+vv2*yy0+vv3*zz0+vv4*ww0+vv5*gg0+vv6*mm0  )+q1*(  yy0      +vv1*xx0+vv2*yy0+vv3*zz0+vv4*ww0+vv5*gg0+vv6*mm0    ))/gamma(q1+2);
yy(1)=yy0+h^q2*((    zz0      +vv7*xx0+vv8*yy0+vv9*zz0+vv10*ww0+vv11*gg0+vv12*mm0   )+q2*(  zz0      +vv7*xx0+vv8*yy0+vv9*zz0+vv10*ww0+vv11*gg0+vv12*mm0     ))/gamma(q2+2);
zz(1)=zz0+h^q3*((    -24*xx0-24*yy0-3*zz0+mm0        +vv13*xx0+vv14*yy0+vv15*zz0+vv16*ww0+vv17*gg0+vv18*mm0   )+q3*(  -24*xx0-24*yy0-3*zz0+mm0        +vv13*xx0+vv14*yy0+vv15*zz0+vv16*ww0+vv17*gg0+vv18*mm0   ))/gamma(q3+2);
ww(1)=ww0+h^q4*((   314*gg0    +vv19*xx0+vv20*yy0+vv21*zz0+vv22*ww0+vv23*gg0+vv24*mm0   )+q4*(  314*gg0    +vv19*xx0+vv20*yy0+vv21*zz0+vv22*ww0+vv23*gg0+vv24*mm0   ))/gamma(q4+2);
gg(1)=gg0+h^q5*((   6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +vv25*xx0+vv26*yy0+vv27*zz0+vv28*ww0+vv29*gg0+vv30*mm0     )+q5*(  6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +vv25*xx0+vv26*yy0+vv27*zz0+vv28*ww0+vv29*gg0+vv30*mm0   ))/gamma(q5+2);
mm(1)=mm0+h^q6*((   -1/Tya*mm0   +vv31*xx0+vv32*yy0+vv33*zz0+vv34*ww0+vv35*gg0+vv36*mm0  )+q6*(  -10*mm0   +vv31*xx0+vv32*yy0+vv33*zz0+vv34*ww0+vv35*gg0+vv36*mm0   ))/gamma(q6+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td2/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(     yy0      +vv1*xx0+vv2*yy0+vv3*zz0+vv4*ww0+vv5*gg0+vv6*mm0     );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     zz0      +vv7*xx0+vv8*yy0+vv9*zz0+vv10*ww0+vv11*gg0+vv12*mm0    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -24*xx0-24*yy0-3*zz0+mm0        +vv13*xx0+vv14*yy0+vv15*zz0+vv16*ww0+vv17*gg0+vv18*mm0  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*gg0    +vv19*xx0+vv20*yy0+vv21*zz0+vv22*ww0+vv23*gg0+vv24*mm0  );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(    6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +vv25*xx0+vv26*yy0+vv27*zz0+vv28*ww0+vv29*gg0+vv30*mm0   );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(   -1/Tya*mm0   +vv31*xx0+vv32*yy0+vv33*zz0+vv34*ww0+vv35*gg0+vv36*mm0    );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(    yy0      +vv1*xx0+vv2*yy0+vv3*zz0+vv4*ww0+vv5*gg0+vv6*mm0  );
N2=((n+1)^q2-n^q2)*(     zz0      +vv7*xx0+vv8*yy0+vv9*zz0+vv10*ww0+vv11*gg0+vv12*mm0     );
N3=((n+1)^q3-n^q3)*(    -24*xx0-24*yy0-3*zz0+mm0        +vv13*xx0+vv14*yy0+vv15*zz0+vv16*ww0+vv17*gg0+vv18*mm0 );
N4=((n+1)^q4-n^q4)*(      314*gg0    +vv19*xx0+vv20*yy0+vv21*zz0+vv22*ww0+vv23*gg0+vv24*mm0     );
N5=((n+1)^q5-n^q5)*(   6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +vv25*xx0+vv26*yy0+vv27*zz0+vv28*ww0+vv29*gg0+vv30*mm0   );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mm0   +vv31*xx0+vv32*yy0+vv33*zz0+vv34*ww0+vv35*gg0+vv36*mm0       );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    yy(j)      +vv1*xx(j)+vv2*yy(j)+vv3*zz(j)+vv4*ww(j)+vv5*gg(j)+vv6*mm(j)   );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     zz(j)      +vv7*xx(j)+vv8*yy(j)+vv9*zz(j)+vv10*ww(j)+vv11*gg(j)+vv12*mm(j)     );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +vv13*xx(j)+vv14*yy(j)+vv15*zz(j)+vv16*ww(j)+vv17*gg(j)+vv18*mm(j)   );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     314*gg(j)    +vv19*xx(j)+vv20*yy(j)+vv21*zz(j)+vv22*ww(j)+vv23*gg(j)+vv24*mm(j)   );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(    6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +vv25*xx(j)+vv26*yy(j)+vv27*zz(j)+vv28*ww(j)+vv29*gg(j)+vv30*mm(j)      );
        M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(  -1/Tya*mm(j)   +vv31*xx(j)+vv32*yy(j)+vv33*zz(j)+vv34*ww(j)+vv35*gg(j)+vv36*mm(j)    );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(      yy(j)      +vv1*xx(j)+vv2*yy(j)+vv3*zz(j)+vv4*ww(j)+vv5*gg(j)+vv6*mm(j)    );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(    zz(j)      +vv7*xx(j)+vv8*yy(j)+vv9*zz(j)+vv10*ww(j)+vv11*gg(j)+vv12*mm(j)  );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(    -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +vv13*xx(j)+vv14*yy(j)+vv15*zz(j)+vv16*ww(j)+vv17*gg(j)+vv18*mm(j)    );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(     314*gg(j)    +vv19*xx(j)+vv20*yy(j)+vv21*zz(j)+vv22*ww(j)+vv23*gg(j)+vv24*mm(j)    );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(    6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +vv25*xx(j)+vv26*yy(j)+vv27*zz(j)+vv28*ww(j)+vv29*gg(j)+vv30*mm(j)     );
        N6=N6+((n-j+1)^q6-(n-j)^q6)*(    -1/Tya*mm(j)   +vv31*xx(j)+vv32*yy(j)+vv33*zz(j)+vv34*ww(j)+vv35*gg(j)+vv36*mm(j)      );
end   
xx1(n+1)=xx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
yy1(n+1)=yy0+h^q2*N2/(gamma(q2)*q2);
zz1(n+1)=zz0+h^q3*N3/(gamma(q3)*q3);
ww1(n+1)=ww0+h^q4*N4/(gamma(q4)*q4);
gg1(n+1)=gg0+h^q5*N5/(gamma(q5)*q5);
mm1(n+1)=mm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
xx(n+2)=xx0+h^q1*(   yy(n+1)      +vv1*xx(n+1)+vv2*yy(n+1)+vv3*zz(n+1)+vv4*ww(n+1)+vv5*gg(n+1)+vv6*mm(n+1)   +M1)/gamma(q1+2);
yy(n+2)=yy0+h^q2*(    zz(n+1)      +vv7*xx(n+1)+vv8*yy(n+1)+vv9*zz(n+1)+vv10*ww(n+1)+vv11*gg(n+1)+vv12*mm(n+1)        +M2)/gamma(q2+2);
zz(n+2)=zz0+h^q3*(     -24*xx(n+1)-24*yy(n+1)-3*zz(n+1)+mm(n+1)        +vv13*xx(n+1)+vv14*yy(n+1)+vv15*zz(n+1)+vv16*ww(n+1)+vv17*gg(n+1)+vv18*mm(n+1)    +M3)/gamma(q3+2);
ww(n+2)=ww0+h^q4*(     314*gg(n+1)    +vv19*xx(n+1)+vv20*yy(n+1)+vv21*zz(n+1)+vv22*ww(n+1)+vv23*gg(n+1)+vv24*mm(n+1)        +M4)/gamma(q4+2);
gg(n+2)=gg0+h^q5*(     6.4*xx(n+1)+0.8*zz(n+1)-3/23*sin(ww(n+1))+320/30429*sin(2*ww(n+1))-2/9*gg(n+1)-7/45*mm(n+1)  +vv25*xx(n+1)+vv26*yy(n+1)+vv27*zz(n+1)+vv28*ww(n+1)+vv29*gg(n+1)+vv30*mm(n+1)    +M5)/gamma(q5+2);
mm(n+2)=mm0+h^q6*(     -1/Tya*mm(n+1)   +vv31*xx(n+1)+vv32*yy(n+1)+vv33*zz(n+1)+vv34*ww(n+1)+vv35*gg(n+1)+vv36*mm(n+1)      +M6)/gamma(q6+2);

end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td2/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(    yy0      +vv1*xx0+vv2*yy0+vv3*zz0+vv4*ww0+vv5*gg0+vv6*mm0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(      zz0      +vv7*xx0+vv8*yy0+vv9*zz0+vv10*ww0+vv11*gg0+vv12*mm0     );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -24*xx0-24*yy0-3*zz0+mm0        +vv13*xx0+vv14*yy0+vv15*zz0+vv16*ww0+vv17*gg0+vv18*mm0   );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*gg0    +vv19*xx0+vv20*yy0+vv21*zz0+vv22*ww0+vv23*gg0+vv24*mm0   );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(     6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +vv25*xx0+vv26*yy0+vv27*zz0+vv28*ww0+vv29*gg0+vv30*mm0  );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(       -1/Tya*mm0   +vv31*xx0+vv32*yy0+vv33*zz0+vv34*ww0+vv35*gg0+vv36*mm0    );

N1=((n+1)^q1-n^q1)*(     yy0      +vv1*xx0+vv2*yy0+vv3*zz0+vv4*ww0+vv5*gg0+vv6*mm0       );
N2=((n+1)^q2-n^q2)*(      zz0      +vv7*xx0+vv8*yy0+vv9*zz0+vv10*ww0+vv11*gg0+vv12*mm0      );
N3=((n+1)^q3-n^q3)*(     -24*xx0-24*yy0-3*zz0+mm0        +vv13*xx0+vv14*yy0+vv15*zz0+vv16*ww0+vv17*gg0+vv18*mm0  );
N4=((n+1)^q4-n^q4)*(    314*gg0    +vv19*xx0+vv20*yy0+vv21*zz0+vv22*ww0+vv23*gg0+vv24*mm0   );
N5=((n+1)^q5-n^q5)*(     6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +vv25*xx0+vv26*yy0+vv27*zz0+vv28*ww0+vv29*gg0+vv30*mm0  );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mm0   +vv31*xx0+vv32*yy0+vv33*zz0+vv34*ww0+vv35*gg0+vv36*mm0  );
for j=1:td2/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace) 
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(      yy(j)      +vv1*xx(j)+vv2*yy(j)+vv3*zz(j)+vv4*ww(j)+vv5*gg(j)+vv6*mm(j)      );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zz(j)      +vv7*xx(j)+vv8*yy(j)+vv9*zz(j)+vv10*ww(j)+vv11*gg(j)+vv12*mm(j)  );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +vv13*xx(j)+vv14*yy(j)+vv15*zz(j)+vv16*ww(j)+vv17*gg(j)+vv18*mm(j)  );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(       314*gg(j)    +vv19*xx(j)+vv20*yy(j)+vv21*zz(j)+vv22*ww(j)+vv23*gg(j)+vv24*mm(j)     );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(     6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +vv25*xx(j)+vv26*yy(j)+vv27*zz(j)+vv28*ww(j)+vv29*gg(j)+vv30*mm(j)  );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(       -1/Tya*mm(j)   +vv31*xx(j)+vv32*yy(j)+vv33*zz(j)+vv34*ww(j)+vv35*gg(j)+vv36*mm(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yy(j)      +vv1*xx(j)+vv2*yy(j)+vv3*zz(j)+vv4*ww(j)+vv5*gg(j)+vv6*mm(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zz(j)      +vv7*xx(j)+vv8*yy(j)+vv9*zz(j)+vv10*ww(j)+vv11*gg(j)+vv12*mm(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +vv13*xx(j)+vv14*yy(j)+vv15*zz(j)+vv16*ww(j)+vv17*gg(j)+vv18*mm(j)  );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(       314*gg(j)    +vv19*xx(j)+vv20*yy(j)+vv21*zz(j)+vv22*ww(j)+vv23*gg(j)+vv24*mm(j)       );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(       6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +vv25*xx(j)+vv26*yy(j)+vv27*zz(j)+vv28*ww(j)+vv29*gg(j)+vv30*mm(j) );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(        -1/Tya*mm(j)   +vv31*xx(j)+vv32*yy(j)+vv33*zz(j)+vv34*ww(j)+vv35*gg(j)+vv36*mm(j)       );
end   
for j=td2/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )       
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     yy(j)      +vv1*xx(j)+vv2*yy(j)+vv3*zz(j)+vv4*ww(j)+vv5*gg(j)+vv6*mm(j)     );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zz(j)      +vv7*xx(j)+vv8*yy(j)+vv9*zz(j)+vv10*ww(j)+vv11*gg(j)+vv12*mm(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xx(j)-24*yy(j)-3*zz(j)+mm(j-td2/h)        +vv13*xx(j)+vv14*yy(j)+vv15*zz(j)+vv16*ww(j)+vv17*gg(j)+vv18*mm(j) );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(      314*gg(j)    +vv19*xx(j)+vv20*yy(j)+vv21*zz(j)+vv22*ww(j)+vv23*gg(j)+vv24*mm(j)      );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(       6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j-td2/h)  +vv25*xx(j)+vv26*yy(j)+vv27*zz(j)+vv28*ww(j)+vv29*gg(j)+vv30*mm(j)   );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(       -1/Tya*mm(j-td2/h)   +vv31*xx(j)+vv32*yy(j)+vv33*zz(j)+vv34*ww(j)+vv35*gg(j)+vv36*mm(j)    );      
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yy(j)      +vv1*xx(j)+vv2*yy(j)+vv3*zz(j)+vv4*ww(j)+vv5*gg(j)+vv6*mm(j)      );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zz(j)      +vv7*xx(j)+vv8*yy(j)+vv9*zz(j)+vv10*ww(j)+vv11*gg(j)+vv12*mm(j)   );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xx(j)-24*yy(j)-3*zz(j)+mm(j-td2/h)        +vv13*xx(j)+vv14*yy(j)+vv15*zz(j)+vv16*ww(j)+vv17*gg(j)+vv18*mm(j)   );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(      314*gg(j)    +vv19*xx(j)+vv20*yy(j)+vv21*zz(j)+vv22*ww(j)+vv23*gg(j)+vv24*mm(j)   );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(      6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j-td2/h)  +vv25*xx(j)+vv26*yy(j)+vv27*zz(j)+vv28*ww(j)+vv29*gg(j)+vv30*mm(j)  );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(       -1/Tya*mm(j-td2/h)   +vv31*xx(j)+vv32*yy(j)+vv33*zz(j)+vv34*ww(j)+vv35*gg(j)+vv36*mm(j)    );      
end
xx1(n+1)=xx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
yy1(n+1)=yy0+h^q2*N2/(gamma(q2)*q2);
zz1(n+1)=zz0+h^q3*N3/(gamma(q2)*q3);
ww1(n+1)=ww0+h^q4*N4/(gamma(q2)*q4);
gg1(n+1)=gg0+h^q5*N5/(gamma(q5)*q5);
mm1(n+1)=mm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
xx(n+1)=xx0+h^q1*(     yy(n+1)      +vv1*xx(n+1)+vv2*yy(n+1)+vv3*zz(n+1)+vv4*ww(n+1)+vv5*gg(n+1)+vv6*mm(n+1)   +M1)/gamma(q1+2);
yy(n+1)=yy0+h^q2*(     zz(n+1)      +vv7*xx(n+1)+vv8*yy(n+1)+vv9*zz(n+1)+vv10*ww(n+1)+vv11*gg(n+1)+vv12*mm(n+1)     +M2)/gamma(q2+2);
zz(n+1)=zz0+h^q3*(     -24*xx(n+1)-24*yy(n+1)-3*zz(n+1)+mm(n+1-td2/h)        +vv13*xx(n+1)+vv14*yy(n+1)+vv15*zz(n+1)+vv16*ww(n+1)+vv17*gg(n+1)+vv18*mm(n+1)    +M3)/gamma(q3+2);
ww(n+1)=ww0+h^q4*(     314*gg(n+1)    +vv19*xx(n+1)+vv20*yy(n+1)+vv21*zz(n+1)+vv22*ww(n+1)+vv23*gg(n+1)+vv24*mm(n+1)      +M4)/gamma(q4+2);
gg(n+1)=gg0+h^q5*(     6.4*xx(n+1)+0.8*zz(n+1)-3/23*sin(ww(n+1))+320/30429*sin(2*ww(n+1))-2/9*gg(n+1)-7/45*mm(n+1-td2/h)  +vv25*xx(n+1)+vv26*yy(n+1)+vv27*zz(n+1)+vv28*ww(n+1)+vv29*gg(n+1)+vv30*mm(n+1)    +M5)/gamma(q5+2);
mm(n+1)=mm0+h^q6*(     -1/Tya*mm(n+1-td2/h)   +vv31*xx(n+1)+vv32*yy(n+1)+vv33*zz(n+1)+vv34*ww(n+1)+vv35*gg(n+1)+vv36*mm(n+1)      +M6)/gamma(q6+2);
end

zhengtii3=diag(0.6*(muk1)+0.4*(muk2));
zhengti3=diag(zhengtii3);
% zhengti3=0.6*(kkk1)+0.4*(kkk2);
vvv1=zhengti3(1,1); vvv2=zhengti3(1,2); vvv3=zhengti3(1,3); vvv4=zhengti3(1,4);vvv5=zhengti3(1,5);vvv6=zhengti3(1,6);
vvv7=zhengti3(2,1); vvv8=zhengti3(2,2); vvv9=zhengti3(2,3); vvv10=zhengti3(2,4);vvv11=zhengti3(2,5);vvv12=zhengti3(2,6);
vvv13=zhengti3(3,1);vvv14=zhengti3(3,2);vvv15=zhengti3(3,3);vvv16=zhengti3(3,4);vvv17=zhengti3(3,5);vvv18=zhengti3(3,6);
vvv19=zhengti3(4,1);vvv20=zhengti3(4,2);vvv21=zhengti3(4,3);vvv22=zhengti3(4,4);vvv23=zhengti3(4,5);vvv24=zhengti3(4,6);
vvv25=zhengti3(5,1);vvv26=zhengti3(5,2);vvv27=zhengti3(5,3);vvv28=zhengti3(5,4);vvv29=zhengti3(5,5);vvv30=zhengti3(5,6);
vvv31=zhengti3(6,1);vvv32=zhengti3(6,2);vvv33=zhengti3(6,3);vvv34=zhengti3(6,4);vvv35=zhengti3(6,5);vvv36=zhengti3(6,6);

q1=0.9;q2=0.9;q3=0.9;q4=0.9;q5=0.9;q6=0.9;  
% q1=0.8;q2=0.8;q3=0.8;q4=0.8;q5=0.8;q6=0.8; 

Tya=0.1; td3=0.1;aw1=0.00;rw1=0.00;fp1=2;


xxx0=0.01;yyy0=0.01;zzz0=0.01;www0=0.01;ggg0=0.01;mmm0=0.01;   %initial value

xxx(N+1)=[0];yyy(N+1)=[0];zzz(N+1)=[0];www(N+1)=[0];ggg(N+1)=[0];mmm(N+1)=[0];  %efficiency need improve
xxx1(N+1)=[0];yyy1(N+1)=[0];zzz1(N+1)=[0];www1(N+1)=[0];ggg1(N+1)=[0];mmm1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
xxx1(1)=xxx0+h^q1*(  yyy0      +vvv1*xxx0+vvv2*yyy0+vvv3*zzz0+vvv4*www0+vvv5*ggg0+vvv6*mmm0   )/(gamma(q1)*q1);
yyy1(1)=yyy0+h^q2*(  zzz0      +vvv7*xxx0+vvv8*yyy0+vvv9*zzz0+vvv10*www0+vvv11*ggg0+vvv12*mmm0   )/(gamma(q2)*q2);
zzz1(1)=zzz0+h^q3*(  -24*xxx0-24*yyy0-3*zzz0+mmm0        +vvv13*xxx0+vvv14*yyy0+vvv15*zzz0+vvv16*www0+vvv17*ggg0+vvv18*mmm0    )/(gamma(q3)*q3);
www1(1)=www0+h^q4*(  314*ggg0    +vvv19*xxx0+vvv20*yyy0+vvv21*zzz0+vvv22*www0+vvv23*ggg0+vvv24*mmm0   )/(gamma(q4)*q4);
ggg1(1)=ggg0+h^q5*(  6.4*xxx0+0.8*zzz0-3/23*sin(www0)+320/30429*sin(2*www0)-2/9*ggg0-7/45*mmm0    +vvv25*xxx0+vvv26*yyy0+vvv27*zzz0+vvv28*www0+vvv29*ggg0+vvv30*mmm0   )/(gamma(q5)*q5);
mmm1(1)=mmm0+h^q6*(  -1/Tya*mmm0   +vvv31*xxx0+vvv32*yyy0+vvv33*zzz0+vvv34*www0+vvv35*ggg0+vvv36*mmm0   )/(gamma(q6)*q6);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
xxx(1)=xxx0+h^q1*((  yyy0      +vvv1*xxx0+vvv2*yyy0+vvv3*zzz0+vvv4*www0+vvv5*ggg0+vvv6*mmm0  )+q1*(  yyy0      +vvv1*xxx0+vvv2*yyy0+vvv3*zzz0+vvv4*www0+vvv5*ggg0+vvv6*mmm0    ))/gamma(q1+2);
yyy(1)=yyy0+h^q2*((    zzz0      +vvv7*xxx0+vvv8*yyy0+vvv9*zzz0+vvv10*www0+vvv11*ggg0+vvv12*mmm0   )+q2*(  zzz0      +vvv7*xxx0+vvv8*yyy0+vvv9*zzz0+vvv10*www0+vvv11*ggg0+vvv12*mmm0     ))/gamma(q2+2);
zzz(1)=zzz0+h^q3*((    -24*xxx0-24*yyy0-3*zzz0+mmm0        +vvv13*xxx0+vvv14*yyy0+vvv15*zzz0+vvv16*www0+vvv17*ggg0+vvv18*mmm0   )+q3*(  -24*xxx0-24*yyy0-3*zzz0+mmm0        +vvv13*xxx0+vvv14*yyy0+vvv15*zzz0+vvv16*www0+vvv17*ggg0+vvv18*mmm0   ))/gamma(q3+2);
www(1)=www0+h^q4*((   314*ggg0    +vvv19*xxx0+vvv20*yyy0+vvv21*zzz0+vvv22*www0+vvv23*ggg0+vvv24*mmm0   )+q4*(  314*ggg0    +vvv19*xxx0+vvv20*yyy0+vvv21*zzz0+vvv22*www0+vvv23*ggg0+vvv24*mmm0   ))/gamma(q4+2);
ggg(1)=ggg0+h^q5*((   6.4*xxx0+0.8*zzz0-3/23*sin(www0)+320/30429*sin(2*www0)-2/9*ggg0-7/45*mmm0  +vvv25*xxx0+vvv26*yyy0+vvv27*zzz0+vvv28*www0+vvv29*ggg0+vvv30*mmm0     )+q5*(  6.4*xxx0+0.8*zzz0-3/23*sin(www0)+320/30429*sin(2*www0)-2/9*ggg0-7/45*mmm0  +vvv25*xxx0+vvv26*yyy0+vvv27*zzz0+vvv28*www0+vvv29*ggg0+vvv30*mmm0   ))/gamma(q5+2);
mmm(1)=mmm0+h^q6*((   -1/Tya*mmm0   +vvv31*xxx0+vvv32*yyy0+vvv33*zzz0+vvv34*www0+vvv35*ggg0+vvv36*mmm0  )+q6*(  -10*mmm0   +vvv31*xxx0+vvv32*yyy0+vvv33*zzz0+vvv34*www0+vvv35*ggg0+vvv36*mmm0   ))/gamma(q6+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td3/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(     yyy0      +vvv1*xxx0+vvv2*yyy0+vvv3*zzz0+vvv4*www0+vvv5*ggg0+vvv6*mmm0     );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     zzz0      +vvv7*xxx0+vvv8*yyy0+vvv9*zzz0+vvv10*www0+vvv11*ggg0+vvv12*mmm0    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -24*xxx0-24*yyy0-3*zzz0+mmm0        +vvv13*xxx0+vvv14*yyy0+vvv15*zzz0+vvv16*www0+vvv17*ggg0+vvv18*mmm0  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*ggg0    +vvv19*xxx0+vvv20*yyy0+vvv21*zzz0+vvv22*www0+vvv23*ggg0+vvv24*mmm0  );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(    6.4*xxx0+0.8*zzz0-3/23*sin(www0)+320/30429*sin(2*www0)-2/9*ggg0-7/45*mmm0  +vvv25*xxx0+vvv26*yyy0+vvv27*zzz0+vvv28*www0+vvv29*ggg0+vvv30*mmm0   );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(   -1/Tya*mmm0   +vvv31*xxx0+vvv32*yyy0+vvv33*zzz0+vvv34*www0+vvv35*ggg0+vvv36*mmm0    );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(    yyy0      +vvv1*xxx0+vvv2*yyy0+vvv3*zzz0+vvv4*www0+vvv5*ggg0+vvv6*mmm0  );
N2=((n+1)^q2-n^q2)*(     zzz0      +vvv7*xxx0+vvv8*yyy0+vvv9*zzz0+vvv10*www0+vvv11*ggg0+vvv12*mmm0     );
N3=((n+1)^q3-n^q3)*(    -24*xxx0-24*yyy0-3*zzz0+mmm0        +vvv13*xxx0+vvv14*yyy0+vvv15*zzz0+vvv16*www0+vvv17*ggg0+vvv18*mmm0 );
N4=((n+1)^q4-n^q4)*(      314*ggg0    +vvv19*xxx0+vvv20*yyy0+vvv21*zzz0+vvv22*www0+vvv23*ggg0+vvv24*mmm0     );
N5=((n+1)^q5-n^q5)*(   6.4*xxx0+0.8*zzz0-3/23*sin(www0)+320/30429*sin(2*www0)-2/9*ggg0-7/45*mmm0  +vvv25*xxx0+vvv26*yyy0+vvv27*zzz0+vvv28*www0+vvv29*ggg0+vvv30*mmm0   );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mmm0   +vvv31*xxx0+vvv32*yyy0+vvv33*zzz0+vvv34*www0+vvv35*ggg0+vvv36*mmm0       );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    yyy(j)      +vvv1*xxx(j)+vvv2*yyy(j)+vvv3*zzz(j)+vvv4*www(j)+vvv5*ggg(j)+vvv6*mmm(j)   );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     zzz(j)      +vvv7*xxx(j)+vvv8*yyy(j)+vvv9*zzz(j)+vvv10*www(j)+vvv11*ggg(j)+vvv12*mmm(j)     );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xxx(j)-24*yyy(j)-3*zzz(j)+mmm(j)        +vvv13*xxx(j)+vvv14*yyy(j)+vvv15*zzz(j)+vvv16*www(j)+vvv17*ggg(j)+vvv18*mmm(j)   );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw1*sin(fp1*j)+rw1*rand(1)+    314*ggg(j)    +vvv19*xxx(j)+vvv20*yyy(j)+vvv21*zzz(j)+vvv22*www(j)+vvv23*ggg(j)+vvv24*mmm(j)   );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(  aw1*sin(fp1*j)+rw1*rand(1)+  6.4*xxx(j)+0.8*zzz(j)-3/23*sin(www(j))+320/30429*sin(2*www(j))-2/9*ggg(j)-7/45*mmm(j)  +vvv25*xxx(j)+vvv26*yyy(j)+vvv27*zzz(j)+vvv28*www(j)+vvv29*ggg(j)+vvv30*mmm(j)      );
        M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*( aw1*sin(fp1*j)+rw1*rand(1)+ -1/Tya*mmm(j)   +vvv31*xxx(j)+vvv32*yyy(j)+vvv33*zzz(j)+vvv34*www(j)+vvv35*ggg(j)+vvv36*mmm(j)    );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(      yyy(j)      +vvv1*xxx(j)+vvv2*yyy(j)+vvv3*zzz(j)+vvv4*www(j)+vvv5*ggg(j)+vvv6*mmm(j)    );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(    zzz(j)      +vvv7*xxx(j)+vvv8*yyy(j)+vvv9*zzz(j)+vvv10*www(j)+vvv11*ggg(j)+vvv12*mmm(j)  );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(    -24*xxx(j)-24*yyy(j)-3*zzz(j)+mmm(j)        +vvv13*xxx(j)+vvv14*yyy(j)+vvv15*zzz(j)+vvv16*www(j)+vvv17*ggg(j)+vvv18*mmm(j)    );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(  aw1*sin(fp1*j)+rw1*rand(1)+   314*ggg(j)    +vvv19*xxx(j)+vvv20*yyy(j)+vvv21*zzz(j)+vvv22*www(j)+vvv23*ggg(j)+vvv24*mmm(j)    );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(  aw1*sin(fp1*j)+rw1*rand(1)+  6.4*xxx(j)+0.8*zzz(j)-3/23*sin(www(j))+320/30429*sin(2*www(j))-2/9*ggg(j)-7/45*mmm(j)  +vvv25*xxx(j)+vvv26*yyy(j)+vvv27*zzz(j)+vvv28*www(j)+vvv29*ggg(j)+vvv30*mmm(j)     );
        N6=N6+((n-j+1)^q6-(n-j)^q6)*(  aw1*sin(fp1*j)+rw1*rand(1)+  -1/Tya*mmm(j)   +vvv31*xxx(j)+vvv32*yyy(j)+vvv33*zzz(j)+vvv34*www(j)+vvv35*ggg(j)+vvv36*mmm(j)      );
end   
xxx1(n+1)=xxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
yyy1(n+1)=yyy0+h^q2*N2/(gamma(q2)*q2);
zzz1(n+1)=zzz0+h^q3*N3/(gamma(q3)*q3);
www1(n+1)=www0+h^q4*N4/(gamma(q4)*q4);
ggg1(n+1)=ggg0+h^q5*N5/(gamma(q5)*q5);
mmm1(n+1)=mmm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
xxx(n+2)=xxx0+h^q1*(   yyy(n+1)      +vvv1*xxx(n+1)+vvv2*yyy(n+1)+vvv3*zzz(n+1)+vvv4*www(n+1)+vvv5*ggg(n+1)+vvv6*mmm(n+1)   +M1)/gamma(q1+2);
yyy(n+2)=yyy0+h^q2*(    zzz(n+1)      +vvv7*xxx(n+1)+vvv8*yyy(n+1)+vvv9*zzz(n+1)+vvv10*www(n+1)+vvv11*ggg(n+1)+vvv12*mmm(n+1)        +M2)/gamma(q2+2);
zzz(n+2)=zzz0+h^q3*(     -24*xxx(n+1)-24*yyy(n+1)-3*zzz(n+1)+mmm(n+1)        +vvv13*xxx(n+1)+vvv14*yyy(n+1)+vvv15*zzz(n+1)+vvv16*www(n+1)+vvv17*ggg(n+1)+vvv18*mmm(n+1)    +M3)/gamma(q3+2);
www(n+2)=www0+h^q4*(     314*ggg(n+1)    +vvv19*xxx(n+1)+vvv20*yyy(n+1)+vvv21*zzz(n+1)+vvv22*www(n+1)+vvv23*ggg(n+1)+vvv24*mmm(n+1)        +M4)/gamma(q4+2);
ggg(n+2)=ggg0+h^q5*(     6.4*xxx(n+1)+0.8*zzz(n+1)-3/23*sin(www(n+1))+320/30429*sin(2*www(n+1))-2/9*ggg(n+1)-7/45*mmm(n+1)  +vvv25*xxx(n+1)+vvv26*yyy(n+1)+vvv27*zzz(n+1)+vvv28*www(n+1)+vvv29*ggg(n+1)+vvv30*mmm(n+1)    +M5)/gamma(q5+2);
mmm(n+2)=mmm0+h^q6*(     -1/Tya*mmm(n+1)   +vvv31*xxx(n+1)+vvv32*yyy(n+1)+vvv33*zzz(n+1)+vvv34*www(n+1)+vvv35*ggg(n+1)+vvv36*mmm(n+1)      +M6)/gamma(q6+2);

end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td3/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(    yyy0      +vvv1*xxx0+vvv2*yyy0+vvv3*zzz0+vvv4*www0+vvv5*ggg0+vvv6*mmm0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(      zzz0      +vvv7*xxx0+vvv8*yyy0+vvv9*zzz0+vvv10*www0+vvv11*ggg0+vvv12*mmm0     );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -24*xxx0-24*yyy0-3*zzz0+mmm0        +vvv13*xxx0+vvv14*yyy0+vvv15*zzz0+vvv16*www0+vvv17*ggg0+vvv18*mmm0   );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*ggg0    +vvv19*xxx0+vvv20*yyy0+vvv21*zzz0+vvv22*www0+vvv23*ggg0+vvv24*mmm0   );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(     6.4*xxx0+0.8*zzz0-3/23*sin(www0)+320/30429*sin(2*www0)-2/9*ggg0-7/45*mmm0  +vvv25*xxx0+vvv26*yyy0+vvv27*zzz0+vvv28*www0+vvv29*ggg0+vvv30*mmm0  );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(       -1/Tya*mmm0   +vvv31*xxx0+vvv32*yyy0+vvv33*zzz0+vvv34*www0+vvv35*ggg0+vvv36*mmm0    );

N1=((n+1)^q1-n^q1)*(     yyy0      +vvv1*xxx0+vvv2*yyy0+vvv3*zzz0+vvv4*www0+vvv5*ggg0+vvv6*mmm0       );
N2=((n+1)^q2-n^q2)*(      zzz0      +vvv7*xxx0+vvv8*yyy0+vvv9*zzz0+vvv10*www0+vvv11*ggg0+vvv12*mmm0      );
N3=((n+1)^q3-n^q3)*(     -24*xxx0-24*yyy0-3*zzz0+mmm0        +vvv13*xxx0+vvv14*yyy0+vvv15*zzz0+vvv16*www0+vvv17*ggg0+vvv18*mmm0  );
N4=((n+1)^q4-n^q4)*(    314*ggg0    +vvv19*xxx0+vvv20*yyy0+vvv21*zzz0+vvv22*www0+vvv23*ggg0+vvv24*mmm0   );
N5=((n+1)^q5-n^q5)*(     6.4*xxx0+0.8*zzz0-3/23*sin(www0)+320/30429*sin(2*www0)-2/9*ggg0-7/45*mmm0  +vvv25*xxx0+vvv26*yyy0+vvv27*zzz0+vvv28*www0+vvv29*ggg0+vvv30*mmm0  );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mmm0   +vvv31*xxx0+vvv32*yyy0+vvv33*zzz0+vvv34*www0+vvv35*ggg0+vvv36*mmm0  );
for j=1:td3/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace) 
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(      yyy(j)      +vvv1*xxx(j)+vvv2*yyy(j)+vvv3*zzz(j)+vvv4*www(j)+vvv5*ggg(j)+vvv6*mmm(j)      );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zzz(j)      +vvv7*xxx(j)+vvv8*yyy(j)+vvv9*zzz(j)+vvv10*www(j)+vvv11*ggg(j)+vvv12*mmm(j)  );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -24*xxx(j)-24*yyy(j)-3*zzz(j)+mmm(j)        +vvv13*xxx(j)+vvv14*yyy(j)+vvv15*zzz(j)+vvv16*www(j)+vvv17*ggg(j)+vvv18*mmm(j)  );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aw1*sin(fp1*j)+rw1*rand(1)+     314*ggg(j)    +vvv19*xxx(j)+vvv20*yyy(j)+vvv21*zzz(j)+vvv22*www(j)+vvv23*ggg(j)+vvv24*mmm(j)     );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(  aw1*sin(fp1*j)+rw1*rand(1)+   6.4*xxx(j)+0.8*zzz(j)-3/23*sin(www(j))+320/30429*sin(2*www(j))-2/9*ggg(j)-7/45*mmm(j)  +vvv25*xxx(j)+vvv26*yyy(j)+vvv27*zzz(j)+vvv28*www(j)+vvv29*ggg(j)+vvv30*mmm(j)  );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(   aw1*sin(fp1*j)+rw1*rand(1)+    -1/Tya*mmm(j)   +vvv31*xxx(j)+vvv32*yyy(j)+vvv33*zzz(j)+vvv34*www(j)+vvv35*ggg(j)+vvv36*mmm(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yyy(j)      +vvv1*xxx(j)+vvv2*yyy(j)+vvv3*zzz(j)+vvv4*www(j)+vvv5*ggg(j)+vvv6*mmm(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zzz(j)      +vvv7*xxx(j)+vvv8*yyy(j)+vvv9*zzz(j)+vvv10*www(j)+vvv11*ggg(j)+vvv12*mmm(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xxx(j)-24*yyy(j)-3*zzz(j)+mmm(j)        +vvv13*xxx(j)+vvv14*yyy(j)+vvv15*zzz(j)+vvv16*www(j)+vvv17*ggg(j)+vvv18*mmm(j)  );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(  aw1*sin(fp1*j)+rw1*rand(1)+     314*ggg(j)    +vvv19*xxx(j)+vvv20*yyy(j)+vvv21*zzz(j)+vvv22*www(j)+vvv23*ggg(j)+vvv24*mmm(j)       );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw1*sin(fp1*j)+rw1*rand(1)+    6.4*xxx(j)+0.8*zzz(j)-3/23*sin(www(j))+320/30429*sin(2*www(j))-2/9*ggg(j)-7/45*mmm(j)  +vvv25*xxx(j)+vvv26*yyy(j)+vvv27*zzz(j)+vvv28*www(j)+vvv29*ggg(j)+vvv30*mmm(j) );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(    aw1*sin(fp1*j)+rw1*rand(1)+    -1/Tya*mmm(j)   +vvv31*xxx(j)+vvv32*yyy(j)+vvv33*zzz(j)+vvv34*www(j)+vvv35*ggg(j)+vvv36*mmm(j)       );
end   
for j=td3/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )       
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     yyy(j)      +vvv1*xxx(j)+vvv2*yyy(j)+vvv3*zzz(j)+vvv4*www(j)+vvv5*ggg(j)+vvv6*mmm(j)     );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zzz(j)      +vvv7*xxx(j)+vvv8*yyy(j)+vvv9*zzz(j)+vvv10*www(j)+vvv11*ggg(j)+vvv12*mmm(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xxx(j)-24*yyy(j)-3*zzz(j)+mmm(j-td3/h)        +vvv13*xxx(j)+vvv14*yyy(j)+vvv15*zzz(j)+vvv16*www(j)+vvv17*ggg(j)+vvv18*mmm(j) );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aw1*sin(fp1*j)+rw1*rand(1)+    314*ggg(j)    +vvv19*xxx(j)+vvv20*yyy(j)+vvv21*zzz(j)+vvv22*www(j)+vvv23*ggg(j)+vvv24*mmm(j)      );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(  aw1*sin(fp1*j)+rw1*rand(1)+     6.4*xxx(j)+0.8*zzz(j)-3/23*sin(www(j))+320/30429*sin(2*www(j))-2/9*ggg(j)-7/45*mmm(j-td3/h)  +vvv25*xxx(j)+vvv26*yyy(j)+vvv27*zzz(j)+vvv28*www(j)+vvv29*ggg(j)+vvv30*mmm(j)   );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(   aw1*sin(fp1*j)+rw1*rand(1)+    -1/Tya*mmm(j-td3/h)   +vvv31*xxx(j)+vvv32*yyy(j)+vvv33*zzz(j)+vvv34*www(j)+vvv35*ggg(j)+vvv36*mmm(j)    );      
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yyy(j)      +vvv1*xxx(j)+vvv2*yyy(j)+vvv3*zzz(j)+vvv4*www(j)+vvv5*ggg(j)+vvv6*mmm(j)      );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zzz(j)      +vvv7*xxx(j)+vvv8*yyy(j)+vvv9*zzz(j)+vvv10*www(j)+vvv11*ggg(j)+vvv12*mmm(j)   );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xxx(j)-24*yyy(j)-3*zzz(j)+mmm(j-td3/h)        +vvv13*xxx(j)+vvv14*yyy(j)+vvv15*zzz(j)+vvv16*www(j)+vvv17*ggg(j)+vvv18*mmm(j)   );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw1*sin(fp1*j)+rw1*rand(1)+   314*ggg(j)    +vvv19*xxx(j)+vvv20*yyy(j)+vvv21*zzz(j)+vvv22*www(j)+vvv23*ggg(j)+vvv24*mmm(j)   );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw1*sin(fp1*j)+rw1*rand(1)+   6.4*xxx(j)+0.8*zzz(j)-3/23*sin(www(j))+320/30429*sin(2*www(j))-2/9*ggg(j)-7/45*mmm(j-td3/h)  +vvv25*xxx(j)+vvv26*yyy(j)+vvv27*zzz(j)+vvv28*www(j)+vvv29*ggg(j)+vvv30*mmm(j)  );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(   aw1*sin(fp1*j)+rw1*rand(1)+    -1/Tya*mmm(j-td3/h)   +vvv31*xxx(j)+vvv32*yyy(j)+vvv33*zzz(j)+vvv34*www(j)+vvv35*ggg(j)+vvv36*mmm(j)    );      
end
xxx1(n+1)=xxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
yyy1(n+1)=yyy0+h^q2*N2/(gamma(q2)*q2);
zzz1(n+1)=zzz0+h^q3*N3/(gamma(q2)*q3);
www1(n+1)=www0+h^q4*N4/(gamma(q2)*q4);
ggg1(n+1)=ggg0+h^q5*N5/(gamma(q5)*q5);
mmm1(n+1)=mmm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
xxx(n+1)=xxx0+h^q1*(     yyy(n+1)      +vvv1*xxx(n+1)+vvv2*yyy(n+1)+vvv3*zzz(n+1)+vvv4*www(n+1)+vvv5*ggg(n+1)+vvv6*mmm(n+1)   +M1)/gamma(q1+2);
yyy(n+1)=yyy0+h^q2*(     zzz(n+1)      +vvv7*xxx(n+1)+vvv8*yyy(n+1)+vvv9*zzz(n+1)+vvv10*www(n+1)+vvv11*ggg(n+1)+vvv12*mmm(n+1)     +M2)/gamma(q2+2);
zzz(n+1)=zzz0+h^q3*(     -24*xxx(n+1)-24*yyy(n+1)-3*zzz(n+1)+mmm(n+1-td3/h)        +vvv13*xxx(n+1)+vvv14*yyy(n+1)+vvv15*zzz(n+1)+vvv16*www(n+1)+vvv17*ggg(n+1)+vvv18*mmm(n+1)    +M3)/gamma(q3+2);
www(n+1)=www0+h^q4*(     314*ggg(n+1)    +vvv19*xxx(n+1)+vvv20*yyy(n+1)+vvv21*zzz(n+1)+vvv22*www(n+1)+vvv23*ggg(n+1)+vvv24*mmm(n+1)      +M4)/gamma(q4+2);
ggg(n+1)=ggg0+h^q5*(     6.4*xxx(n+1)+0.8*zzz(n+1)-3/23*sin(www(n+1))+320/30429*sin(2*www(n+1))-2/9*ggg(n+1)-7/45*mmm(n+1-td3/h)  +vvv25*xxx(n+1)+vvv26*yyy(n+1)+vvv27*zzz(n+1)+vvv28*www(n+1)+vvv29*ggg(n+1)+vvv30*mmm(n+1)    +M5)/gamma(q5+2);
mmm(n+1)=mmm0+h^q6*(     -1/Tya*mmm(n+1-td3/h)   +vvv31*xxx(n+1)+vvv32*yyy(n+1)+vvv33*zzz(n+1)+vvv34*www(n+1)+vvv35*ggg(n+1)+vvv36*mmm(n+1)      +M6)/gamma(q6+2);
end


zhengtii4=diag(0.6*(muk1)+0.4*(muk2));
zhengti4=diag(zhengtii4);
% zhengti4=0.6*(kkkk1)+0.4*(kkkk2);
vvvv1=zhengti4(1,1); vvvv2=zhengti4(1,2); vvvv3=zhengti4(1,3); vvvv4=zhengti4(1,4);vvvv5=zhengti4(1,5);vvvv6=zhengti4(1,6);
vvvv7=zhengti4(2,1); vvvv8=zhengti4(2,2); vvvv9=zhengti4(2,3); vvvv10=zhengti4(2,4);vvvv11=zhengti4(2,5);vvvv12=zhengti4(2,6);
vvvv13=zhengti4(3,1);vvvv14=zhengti4(3,2);vvvv15=zhengti4(3,3);vvvv16=zhengti4(3,4);vvvv17=zhengti4(3,5);vvvv18=zhengti4(3,6);
vvvv19=zhengti4(4,1);vvvv20=zhengti4(4,2);vvvv21=zhengti4(4,3);vvvv22=zhengti4(4,4);vvvv23=zhengti4(4,5);vvvv24=zhengti4(4,6);
vvvv25=zhengti4(5,1);vvvv26=zhengti4(5,2);vvvv27=zhengti4(5,3);vvvv28=zhengti4(5,4);vvvv29=zhengti4(5,5);vvvv30=zhengti4(5,6);
vvvv31=zhengti4(6,1);vvvv32=zhengti4(6,2);vvvv33=zhengti4(6,3);vvvv34=zhengti4(6,4);vvvv35=zhengti4(6,5);vvvv36=zhengti4(6,6);

q1=0.8;q2=0.8;q3=0.8;q4=0.8;q5=0.8;q6=0.8; 
% q1=0.75;q2=0.75;q3=0.75;q4=0.75;q5=0.75;q6=0.75;  

Tya=0.1; td4=0.1;aw2=0.00;rw2=0.00;fp2=1000;


xxxx0=0.01;yyyy0=0.01;zzzz0=0.01;wwww0=0.01;gggg0=0.01;mmmm0=0.01;   %initial value

xxxx(N+1)=[0];yyyy(N+1)=[0];zzzz(N+1)=[0];wwww(N+1)=[0];gggg(N+1)=[0];mmmm(N+1)=[0];  %efficiency need improve
xxxx1(N+1)=[0];yyyy1(N+1)=[0];zzzz1(N+1)=[0];wwww1(N+1)=[0];gggg1(N+1)=[0];mmmm1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
xxxx1(1)=xxxx0+h^q1*(  yyyy0      +vvvv1*xxxx0+vvvv2*yyyy0+vvvv3*zzzz0+vvvv4*wwww0+vvvv5*gggg0+vvvv6*mmmm0   )/(gamma(q1)*q1);
yyyy1(1)=yyyy0+h^q2*(  zzzz0      +vvvv7*xxxx0+vvvv8*yyyy0+vvvv9*zzzz0+vvvv10*wwww0+vvvv11*gggg0+vvvv12*mmmm0   )/(gamma(q2)*q2);
zzzz1(1)=zzzz0+h^q3*(  -24*xxxx0-24*yyyy0-3*zzzz0+mmmm0        +vvvv13*xxxx0+vvvv14*yyyy0+vvvv15*zzzz0+vvvv16*wwww0+vvvv17*gggg0+vvvv18*mmmm0    )/(gamma(q3)*q3);
wwww1(1)=wwww0+h^q4*(  314*gggg0    +vvvv19*xxxx0+vvvv20*yyyy0+vvvv21*zzzz0+vvvv22*wwww0+vvvv23*gggg0+vvvv24*mmmm0   )/(gamma(q4)*q4);
gggg1(1)=gggg0+h^q5*(  6.4*xxxx0+0.8*zzzz0-3/23*sin(wwww0)+320/30429*sin(2*wwww0)-2/9*gggg0-7/45*mmmm0    +vvvv25*xxxx0+vvvv26*yyyy0+vvvv27*zzzz0+vvvv28*wwww0+vvvv29*gggg0+vvvv30*mmmm0   )/(gamma(q5)*q5);
mmmm1(1)=mmmm0+h^q6*(  -1/Tya*mmmm0   +vvvv31*xxxx0+vvvv32*yyyy0+vvvv33*zzzz0+vvvv34*wwww0+vvvv35*gggg0+vvvv36*mmmm0   )/(gamma(q6)*q6);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
xxxx(1)=xxxx0+h^q1*((  yyyy0      +vvvv1*xxxx0+vvvv2*yyyy0+vvvv3*zzzz0+vvvv4*wwww0+vvvv5*gggg0+vvvv6*mmmm0  )+q1*(  yyyy0      +vvvv1*xxxx0+vvvv2*yyyy0+vvvv3*zzzz0+vvvv4*wwww0+vvvv5*gggg0+vvvv6*mmmm0    ))/gamma(q1+2);
yyyy(1)=yyyy0+h^q2*((    zzzz0      +vvvv7*xxxx0+vvvv8*yyyy0+vvvv9*zzzz0+vvvv10*wwww0+vvvv11*gggg0+vvvv12*mmmm0   )+q2*(  zzzz0      +vvvv7*xxxx0+vvvv8*yyyy0+vvvv9*zzzz0+vvvv10*wwww0+vvvv11*gggg0+vvvv12*mmmm0     ))/gamma(q2+2);
zzzz(1)=zzzz0+h^q3*((    -24*xxxx0-24*yyyy0-3*zzzz0+mmmm0        +vvvv13*xxxx0+vvvv14*yyyy0+vvvv15*zzzz0+vvvv16*wwww0+vvvv17*gggg0+vvvv18*mmmm0   )+q3*(  -24*xxxx0-24*yyyy0-3*zzzz0+mmmm0        +vvvv13*xxxx0+vvvv14*yyyy0+vvvv15*zzzz0+vvvv16*wwww0+vvvv17*gggg0+vvvv18*mmmm0   ))/gamma(q3+2);
wwww(1)=wwww0+h^q4*((   314*gggg0    +vvvv19*xxxx0+vvvv20*yyyy0+vvvv21*zzzz0+vvvv22*wwww0+vvvv23*gggg0+vvvv24*mmmm0   )+q4*(  314*gggg0    +vvvv19*xxxx0+vvvv20*yyyy0+vvvv21*zzzz0+vvvv22*wwww0+vvvv23*gggg0+vvvv24*mmmm0   ))/gamma(q4+2);
gggg(1)=gggg0+h^q5*((   6.4*xxxx0+0.8*zzzz0-3/23*sin(wwww0)+320/30429*sin(2*wwww0)-2/9*gggg0-7/45*mmmm0  +vvvv25*xxxx0+vvvv26*yyyy0+vvvv27*zzzz0+vvvv28*wwww0+vvvv29*gggg0+vvvv30*mmmm0     )+q5*(  6.4*xxxx0+0.8*zzzz0-3/23*sin(wwww0)+320/30429*sin(2*wwww0)-2/9*gggg0-7/45*mmmm0  +vvvv25*xxxx0+vvvv26*yyyy0+vvvv27*zzzz0+vvvv28*wwww0+vvvv29*gggg0+vvvv30*mmmm0   ))/gamma(q5+2);
mmmm(1)=mmmm0+h^q6*((   -1/Tya*mmmm0   +vvvv31*xxxx0+vvvv32*yyyy0+vvvv33*zzzz0+vvvv34*wwww0+vvvv35*gggg0+vvvv36*mmmm0  )+q6*(  -10*mmmm0   +vvvv31*xxxx0+vvvv32*yyyy0+vvvv33*zzzz0+vvvv34*wwww0+vvvv35*gggg0+vvvv36*mmmm0   ))/gamma(q6+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td4/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(     yyyy0      +vvvv1*xxxx0+vvvv2*yyyy0+vvvv3*zzzz0+vvvv4*wwww0+vvvv5*gggg0+vvvv6*mmmm0     );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     zzzz0      +vvvv7*xxxx0+vvvv8*yyyy0+vvvv9*zzzz0+vvvv10*wwww0+vvvv11*gggg0+vvvv12*mmmm0    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -24*xxxx0-24*yyyy0-3*zzzz0+mmmm0        +vvvv13*xxxx0+vvvv14*yyyy0+vvvv15*zzzz0+vvvv16*wwww0+vvvv17*gggg0+vvvv18*mmmm0  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*gggg0    +vvvv19*xxxx0+vvvv20*yyyy0+vvvv21*zzzz0+vvvv22*wwww0+vvvv23*gggg0+vvvv24*mmmm0  );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(    6.4*xxxx0+0.8*zzzz0-3/23*sin(wwww0)+320/30429*sin(2*wwww0)-2/9*gggg0-7/45*mmmm0  +vvvv25*xxxx0+vvvv26*yyyy0+vvvv27*zzzz0+vvvv28*wwww0+vvvv29*gggg0+vvvv30*mmmm0   );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(   -1/Tya*mmmm0   +vvvv31*xxxx0+vvvv32*yyyy0+vvvv33*zzzz0+vvvv34*wwww0+vvvv35*gggg0+vvvv36*mmmm0    );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(    yyyy0      +vvvv1*xxxx0+vvvv2*yyyy0+vvvv3*zzzz0+vvvv4*wwww0+vvvv5*gggg0+vvvv6*mmmm0  );
N2=((n+1)^q2-n^q2)*(     zzzz0      +vvvv7*xxxx0+vvvv8*yyyy0+vvvv9*zzzz0+vvvv10*wwww0+vvvv11*gggg0+vvvv12*mmmm0     );
N3=((n+1)^q3-n^q3)*(    -24*xxxx0-24*yyyy0-3*zzzz0+mmmm0        +vvvv13*xxxx0+vvvv14*yyyy0+vvvv15*zzzz0+vvvv16*wwww0+vvvv17*gggg0+vvvv18*mmmm0 );
N4=((n+1)^q4-n^q4)*(      314*gggg0    +vvvv19*xxxx0+vvvv20*yyyy0+vvvv21*zzzz0+vvvv22*wwww0+vvvv23*gggg0+vvvv24*mmmm0     );
N5=((n+1)^q5-n^q5)*(   6.4*xxxx0+0.8*zzzz0-3/23*sin(wwww0)+320/30429*sin(2*wwww0)-2/9*gggg0-7/45*mmmm0  +vvvv25*xxxx0+vvvv26*yyyy0+vvvv27*zzzz0+vvvv28*wwww0+vvvv29*gggg0+vvvv30*mmmm0   );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mmmm0   +vvvv31*xxxx0+vvvv32*yyyy0+vvvv33*zzzz0+vvvv34*wwww0+vvvv35*gggg0+vvvv36*mmmm0       );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    yyyy(j)      +vvvv1*xxxx(j)+vvvv2*yyyy(j)+vvvv3*zzzz(j)+vvvv4*wwww(j)+vvvv5*gggg(j)+vvvv6*mmmm(j)   );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     zzzz(j)      +vvvv7*xxxx(j)+vvvv8*yyyy(j)+vvvv9*zzzz(j)+vvvv10*wwww(j)+vvvv11*gggg(j)+vvvv12*mmmm(j)     );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xxxx(j)-24*yyyy(j)-3*zzzz(j)+mmmm(j)        +vvvv13*xxxx(j)+vvvv14*yyyy(j)+vvvv15*zzzz(j)+vvvv16*wwww(j)+vvvv17*gggg(j)+vvvv18*mmmm(j)   );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+  314*gggg(j)    +vvvv19*xxxx(j)+vvvv20*yyyy(j)+vvvv21*zzzz(j)+vvvv22*wwww(j)+vvvv23*gggg(j)+vvvv24*mmmm(j)   );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(  aw2*sin(fp2*j)+rw2*rand(1)+  6.4*xxxx(j)+0.8*zzzz(j)-3/23*sin(wwww(j))+320/30429*sin(2*wwww(j))-2/9*gggg(j)-7/45*mmmm(j)  +vvvv25*xxxx(j)+vvvv26*yyyy(j)+vvvv27*zzzz(j)+vvvv28*wwww(j)+vvvv29*gggg(j)+vvvv30*mmmm(j)      );
        M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(  aw2*sin(fp2*j)+rw2*rand(1)+-1/Tya*mmmm(j)   +vvvv31*xxxx(j)+vvvv32*yyyy(j)+vvvv33*zzzz(j)+vvvv34*wwww(j)+vvvv35*gggg(j)+vvvv36*mmmm(j)    );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(      yyyy(j)      +vvvv1*xxxx(j)+vvvv2*yyyy(j)+vvvv3*zzzz(j)+vvvv4*wwww(j)+vvvv5*gggg(j)+vvvv6*mmmm(j)    );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(    zzzz(j)      +vvvv7*xxxx(j)+vvvv8*yyyy(j)+vvvv9*zzzz(j)+vvvv10*wwww(j)+vvvv11*gggg(j)+vvvv12*mmmm(j)  );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(    -24*xxxx(j)-24*yyyy(j)-3*zzzz(j)+mmmm(j)        +vvvv13*xxxx(j)+vvvv14*yyyy(j)+vvvv15*zzzz(j)+vvvv16*wwww(j)+vvvv17*gggg(j)+vvvv18*mmmm(j)    );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw2*sin(fp2*j)+rw2*rand(1)+  314*gggg(j)    +vvvv19*xxxx(j)+vvvv20*yyyy(j)+vvvv21*zzzz(j)+vvvv22*wwww(j)+vvvv23*gggg(j)+vvvv24*mmmm(j)    );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw2*sin(fp2*j)+rw2*rand(1)+ 6.4*xxxx(j)+0.8*zzzz(j)-3/23*sin(wwww(j))+320/30429*sin(2*wwww(j))-2/9*gggg(j)-7/45*mmmm(j)  +vvvv25*xxxx(j)+vvvv26*yyyy(j)+vvvv27*zzzz(j)+vvvv28*wwww(j)+vvvv29*gggg(j)+vvvv30*mmmm(j)     );
        N6=N6+((n-j+1)^q6-(n-j)^q6)*(   aw2*sin(fp2*j)+rw2*rand(1)+ -1/Tya*mmmm(j)   +vvvv31*xxxx(j)+vvvv32*yyyy(j)+vvvv33*zzzz(j)+vvvv34*wwww(j)+vvvv35*gggg(j)+vvvv36*mmmm(j)      );
end   
xxxx1(n+1)=xxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
yyyy1(n+1)=yyyy0+h^q2*N2/(gamma(q2)*q2);
zzzz1(n+1)=zzzz0+h^q3*N3/(gamma(q3)*q3);
wwww1(n+1)=wwww0+h^q4*N4/(gamma(q4)*q4);
gggg1(n+1)=gggg0+h^q5*N5/(gamma(q5)*q5);
mmmm1(n+1)=mmmm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
xxxx(n+2)=xxxx0+h^q1*(   yyyy(n+1)      +vvvv1*xxxx(n+1)+vvvv2*yyyy(n+1)+vvvv3*zzzz(n+1)+vvvv4*wwww(n+1)+vvvv5*gggg(n+1)+vvvv6*mmmm(n+1)   +M1)/gamma(q1+2);
yyyy(n+2)=yyyy0+h^q2*(    zzzz(n+1)      +vvvv7*xxxx(n+1)+vvvv8*yyyy(n+1)+vvvv9*zzzz(n+1)+vvvv10*wwww(n+1)+vvvv11*gggg(n+1)+vvvv12*mmmm(n+1)        +M2)/gamma(q2+2);
zzzz(n+2)=zzzz0+h^q3*(     -24*xxxx(n+1)-24*yyyy(n+1)-3*zzzz(n+1)+mmmm(n+1)        +vvvv13*xxxx(n+1)+vvvv14*yyyy(n+1)+vvvv15*zzzz(n+1)+vvvv16*wwww(n+1)+vvvv17*gggg(n+1)+vvvv18*mmmm(n+1)    +M3)/gamma(q3+2);
wwww(n+2)=wwww0+h^q4*(     314*gggg(n+1)    +vvvv19*xxxx(n+1)+vvvv20*yyyy(n+1)+vvvv21*zzzz(n+1)+vvvv22*wwww(n+1)+vvvv23*gggg(n+1)+vvvv24*mmmm(n+1)        +M4)/gamma(q4+2);
gggg(n+2)=gggg0+h^q5*(     6.4*xxxx(n+1)+0.8*zzzz(n+1)-3/23*sin(wwww(n+1))+320/30429*sin(2*wwww(n+1))-2/9*gggg(n+1)-7/45*mmmm(n+1)  +vvvv25*xxxx(n+1)+vvvv26*yyyy(n+1)+vvvv27*zzzz(n+1)+vvvv28*wwww(n+1)+vvvv29*gggg(n+1)+vvvv30*mmmm(n+1)    +M5)/gamma(q5+2);
mmmm(n+2)=mmmm0+h^q6*(     -1/Tya*mmmm(n+1)   +vvvv31*xxxx(n+1)+vvvv32*yyyy(n+1)+vvvv33*zzzz(n+1)+vvvv34*wwww(n+1)+vvvv35*gggg(n+1)+vvvv36*mmmm(n+1)      +M6)/gamma(q6+2);

end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td4/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(    yyyy0      +vvvv1*xxxx0+vvvv2*yyyy0+vvvv3*zzzz0+vvvv4*wwww0+vvvv5*gggg0+vvvv6*mmmm0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(      zzzz0      +vvvv7*xxxx0+vvvv8*yyyy0+vvvv9*zzzz0+vvvv10*wwww0+vvvv11*gggg0+vvvv12*mmmm0     );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -24*xxxx0-24*yyyy0-3*zzzz0+mmmm0        +vvvv13*xxxx0+vvvv14*yyyy0+vvvv15*zzzz0+vvvv16*wwww0+vvvv17*gggg0+vvvv18*mmmm0   );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*gggg0    +vvvv19*xxxx0+vvvv20*yyyy0+vvvv21*zzzz0+vvvv22*wwww0+vvvv23*gggg0+vvvv24*mmmm0   );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(     6.4*xxxx0+0.8*zzzz0-3/23*sin(wwww0)+320/30429*sin(2*wwww0)-2/9*gggg0-7/45*mmmm0  +vvvv25*xxxx0+vvvv26*yyyy0+vvvv27*zzzz0+vvvv28*wwww0+vvvv29*gggg0+vvvv30*mmmm0  );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(       -1/Tya*mmmm0   +vvvv31*xxxx0+vvvv32*yyyy0+vvvv33*zzzz0+vvvv34*wwww0+vvvv35*gggg0+vvvv36*mmmm0    );

N1=((n+1)^q1-n^q1)*(     yyyy0      +vvvv1*xxxx0+vvvv2*yyyy0+vvvv3*zzzz0+vvvv4*wwww0+vvvv5*gggg0+vvvv6*mmmm0       );
N2=((n+1)^q2-n^q2)*(      zzzz0      +vvvv7*xxxx0+vvvv8*yyyy0+vvvv9*zzzz0+vvvv10*wwww0+vvvv11*gggg0+vvvv12*mmmm0      );
N3=((n+1)^q3-n^q3)*(     -24*xxxx0-24*yyyy0-3*zzzz0+mmmm0        +vvvv13*xxxx0+vvvv14*yyyy0+vvvv15*zzzz0+vvvv16*wwww0+vvvv17*gggg0+vvvv18*mmmm0  );
N4=((n+1)^q4-n^q4)*(    314*gggg0    +vvvv19*xxxx0+vvvv20*yyyy0+vvvv21*zzzz0+vvvv22*wwww0+vvvv23*gggg0+vvvv24*mmmm0   );
N5=((n+1)^q5-n^q5)*(     6.4*xxxx0+0.8*zzzz0-3/23*sin(wwww0)+320/30429*sin(2*wwww0)-2/9*gggg0-7/45*mmmm0  +vvvv25*xxxx0+vvvv26*yyyy0+vvvv27*zzzz0+vvvv28*wwww0+vvvv29*gggg0+vvvv30*mmmm0  );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mmmm0   +vvvv31*xxxx0+vvvv32*yyyy0+vvvv33*zzzz0+vvvv34*wwww0+vvvv35*gggg0+vvvv36*mmmm0  );
for j=1:td4/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace) 
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(      yyyy(j)      +vvvv1*xxxx(j)+vvvv2*yyyy(j)+vvvv3*zzzz(j)+vvvv4*wwww(j)+vvvv5*gggg(j)+vvvv6*mmmm(j)      );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zzzz(j)      +vvvv7*xxxx(j)+vvvv8*yyyy(j)+vvvv9*zzzz(j)+vvvv10*wwww(j)+vvvv11*gggg(j)+vvvv12*mmmm(j)  );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -24*xxxx(j)-24*yyyy(j)-3*zzzz(j)+mmmm(j)        +vvvv13*xxxx(j)+vvvv14*yyyy(j)+vvvv15*zzzz(j)+vvvv16*wwww(j)+vvvv17*gggg(j)+vvvv18*mmmm(j)  );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aw2*sin(fp2*j)+rw2*rand(1)+     314*gggg(j)    +vvvv19*xxxx(j)+vvvv20*yyyy(j)+vvvv21*zzzz(j)+vvvv22*wwww(j)+vvvv23*gggg(j)+vvvv24*mmmm(j)     );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+  6.4*xxxx(j)+0.8*zzzz(j)-3/23*sin(wwww(j))+320/30429*sin(2*wwww(j))-2/9*gggg(j)-7/45*mmmm(j)  +vvvv25*xxxx(j)+vvvv26*yyyy(j)+vvvv27*zzzz(j)+vvvv28*wwww(j)+vvvv29*gggg(j)+vvvv30*mmmm(j)  );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+    -1/Tya*mmmm(j)   +vvvv31*xxxx(j)+vvvv32*yyyy(j)+vvvv33*zzzz(j)+vvvv34*wwww(j)+vvvv35*gggg(j)+vvvv36*mmmm(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yyyy(j)      +vvvv1*xxxx(j)+vvvv2*yyyy(j)+vvvv3*zzzz(j)+vvvv4*wwww(j)+vvvv5*gggg(j)+vvvv6*mmmm(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zzzz(j)      +vvvv7*xxxx(j)+vvvv8*yyyy(j)+vvvv9*zzzz(j)+vvvv10*wwww(j)+vvvv11*gggg(j)+vvvv12*mmmm(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xxxx(j)-24*yyyy(j)-3*zzzz(j)+mmmm(j)        +vvvv13*xxxx(j)+vvvv14*yyyy(j)+vvvv15*zzzz(j)+vvvv16*wwww(j)+vvvv17*gggg(j)+vvvv18*mmmm(j)  );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw2*sin(fp2*j)+rw2*rand(1)+    314*gggg(j)    +vvvv19*xxxx(j)+vvvv20*yyyy(j)+vvvv21*zzzz(j)+vvvv22*wwww(j)+vvvv23*gggg(j)+vvvv24*mmmm(j)       );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw2*sin(fp2*j)+rw2*rand(1)+    6.4*xxxx(j)+0.8*zzzz(j)-3/23*sin(wwww(j))+320/30429*sin(2*wwww(j))-2/9*gggg(j)-7/45*mmmm(j)  +vvvv25*xxxx(j)+vvvv26*yyyy(j)+vvvv27*zzzz(j)+vvvv28*wwww(j)+vvvv29*gggg(j)+vvvv30*mmmm(j) );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(   aw2*sin(fp2*j)+rw2*rand(1)+     -1/Tya*mmmm(j)   +vvvv31*xxxx(j)+vvvv32*yyyy(j)+vvvv33*zzzz(j)+vvvv34*wwww(j)+vvvv35*gggg(j)+vvvv36*mmmm(j)       );
end   
for j=td4/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )       
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     yyyy(j)      +vvvv1*xxxx(j)+vvvv2*yyyy(j)+vvvv3*zzzz(j)+vvvv4*wwww(j)+vvvv5*gggg(j)+vvvv6*mmmm(j)     );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zzzz(j)      +vvvv7*xxxx(j)+vvvv8*yyyy(j)+vvvv9*zzzz(j)+vvvv10*wwww(j)+vvvv11*gggg(j)+vvvv12*mmmm(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xxxx(j)-24*yyyy(j)-3*zzzz(j)+mmmm(j-td4/h)        +vvvv13*xxxx(j)+vvvv14*yyyy(j)+vvvv15*zzzz(j)+vvvv16*wwww(j)+vvvv17*gggg(j)+vvvv18*mmmm(j) );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+   314*gggg(j)    +vvvv19*xxxx(j)+vvvv20*yyyy(j)+vvvv21*zzzz(j)+vvvv22*wwww(j)+vvvv23*gggg(j)+vvvv24*mmmm(j)      );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+    6.4*xxxx(j)+0.8*zzzz(j)-3/23*sin(wwww(j))+320/30429*sin(2*wwww(j))-2/9*gggg(j)-7/45*mmmm(j-td4/h)  +vvvv25*xxxx(j)+vvvv26*yyyy(j)+vvvv27*zzzz(j)+vvvv28*wwww(j)+vvvv29*gggg(j)+vvvv30*mmmm(j)   );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(    aw2*sin(fp2*j)+rw2*rand(1)+   -1/Tya*mmmm(j-td4/h)   +vvvv31*xxxx(j)+vvvv32*yyyy(j)+vvvv33*zzzz(j)+vvvv34*wwww(j)+vvvv35*gggg(j)+vvvv36*mmmm(j)    );      
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yyyy(j)      +vvvv1*xxxx(j)+vvvv2*yyyy(j)+vvvv3*zzzz(j)+vvvv4*wwww(j)+vvvv5*gggg(j)+vvvv6*mmmm(j)      );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zzzz(j)      +vvvv7*xxxx(j)+vvvv8*yyyy(j)+vvvv9*zzzz(j)+vvvv10*wwww(j)+vvvv11*gggg(j)+vvvv12*mmmm(j)   );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xxxx(j)-24*yyyy(j)-3*zzzz(j)+mmmm(j-td4/h)        +vvvv13*xxxx(j)+vvvv14*yyyy(j)+vvvv15*zzzz(j)+vvvv16*wwww(j)+vvvv17*gggg(j)+vvvv18*mmmm(j)   );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw2*sin(fp2*j)+rw2*rand(1)+   314*gggg(j)    +vvvv19*xxxx(j)+vvvv20*yyyy(j)+vvvv21*zzzz(j)+vvvv22*wwww(j)+vvvv23*gggg(j)+vvvv24*mmmm(j)   );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw2*sin(fp2*j)+rw2*rand(1)+   6.4*xxxx(j)+0.8*zzzz(j)-3/23*sin(wwww(j))+320/30429*sin(2*wwww(j))-2/9*gggg(j)-7/45*mmmm(j-td4/h)  +vvvv25*xxxx(j)+vvvv26*yyyy(j)+vvvv27*zzzz(j)+vvvv28*wwww(j)+vvvv29*gggg(j)+vvvv30*mmmm(j)  );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(    aw2*sin(fp2*j)+rw2*rand(1)+   -1/Tya*mmmm(j-td4/h)   +vvvv31*xxxx(j)+vvvv32*yyyy(j)+vvvv33*zzzz(j)+vvvv34*wwww(j)+vvvv35*gggg(j)+vvvv36*mmmm(j)    );      
end
xxxx1(n+1)=xxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
yyyy1(n+1)=yyyy0+h^q2*N2/(gamma(q2)*q2);
zzzz1(n+1)=zzzz0+h^q3*N3/(gamma(q2)*q3);
wwww1(n+1)=wwww0+h^q4*N4/(gamma(q2)*q4);
gggg1(n+1)=gggg0+h^q5*N5/(gamma(q5)*q5);
mmmm1(n+1)=mmmm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
xxxx(n+1)=xxxx0+h^q1*(     yyyy(n+1)      +vvvv1*xxxx(n+1)+vvvv2*yyyy(n+1)+vvvv3*zzzz(n+1)+vvvv4*wwww(n+1)+vvvv5*gggg(n+1)+vvvv6*mmmm(n+1)   +M1)/gamma(q1+2);
yyyy(n+1)=yyyy0+h^q2*(     zzzz(n+1)      +vvvv7*xxxx(n+1)+vvvv8*yyyy(n+1)+vvvv9*zzzz(n+1)+vvvv10*wwww(n+1)+vvvv11*gggg(n+1)+vvvv12*mmmm(n+1)     +M2)/gamma(q2+2);
zzzz(n+1)=zzzz0+h^q3*(     -24*xxxx(n+1)-24*yyyy(n+1)-3*zzzz(n+1)+mmmm(n+1-td4/h)        +vvvv13*xxxx(n+1)+vvvv14*yyyy(n+1)+vvvv15*zzzz(n+1)+vvvv16*wwww(n+1)+vvvv17*gggg(n+1)+vvvv18*mmmm(n+1)    +M3)/gamma(q3+2);
wwww(n+1)=wwww0+h^q4*(     314*gggg(n+1)    +vvvv19*xxxx(n+1)+vvvv20*yyyy(n+1)+vvvv21*zzzz(n+1)+vvvv22*wwww(n+1)+vvvv23*gggg(n+1)+vvvv24*mmmm(n+1)      +M4)/gamma(q4+2);
gggg(n+1)=gggg0+h^q5*(     6.4*xxxx(n+1)+0.8*zzzz(n+1)-3/23*sin(wwww(n+1))+320/30429*sin(2*wwww(n+1))-2/9*gggg(n+1)-7/45*mmmm(n+1-td4/h)  +vvvv25*xxxx(n+1)+vvvv26*yyyy(n+1)+vvvv27*zzzz(n+1)+vvvv28*wwww(n+1)+vvvv29*gggg(n+1)+vvvv30*mmmm(n+1)    +M5)/gamma(q5+2);
mmmm(n+1)=mmmm0+h^q6*(     -1/Tya*mmmm(n+1-td4/h)   +vvvv31*xxxx(n+1)+vvvv32*yyyy(n+1)+vvvv33*zzzz(n+1)+vvvv34*wwww(n+1)+vvvv35*gggg(n+1)+vvvv36*mmmm(n+1)      +M6)/gamma(q6+2);
end

zhengtii4=diag(0.6*(muk1)+0.4*(muk2));
zhengti4=diag(zhengtii4);
% zhengti4=0.6*(kkkk1)+0.4*(kkkk2);
vvvv1=zhengti4(1,1); vvvv2=zhengti4(1,2); vvvv3=zhengti4(1,3); vvvv4=zhengti4(1,4);vvvv5=zhengti4(1,5);vvvv6=zhengti4(1,6);
vvvv7=zhengti4(2,1); vvvv8=zhengti4(2,2); vvvv9=zhengti4(2,3); vvvv10=zhengti4(2,4);vvvv11=zhengti4(2,5);vvvv12=zhengti4(2,6);
vvvv13=zhengti4(3,1);vvvv14=zhengti4(3,2);vvvv15=zhengti4(3,3);vvvv16=zhengti4(3,4);vvvv17=zhengti4(3,5);vvvv18=zhengti4(3,6);
vvvv19=zhengti4(4,1);vvvv20=zhengti4(4,2);vvvv21=zhengti4(4,3);vvvv22=zhengti4(4,4);vvvv23=zhengti4(4,5);vvvv24=zhengti4(4,6);
vvvv25=zhengti4(5,1);vvvv26=zhengti4(5,2);vvvv27=zhengti4(5,3);vvvv28=zhengti4(5,4);vvvv29=zhengti4(5,5);vvvv30=zhengti4(5,6);
vvvv31=zhengti4(6,1);vvvv32=zhengti4(6,2);vvvv33=zhengti4(6,3);vvvv34=zhengti4(6,4);vvvv35=zhengti4(6,5);vvvv36=zhengti4(6,6);

% q1=0.78;q2=0.78;q3=0.78;q4=0.78;q5=0.78;q6=0.78; 
q1=0.74;q2=0.74;q3=0.74;q4=0.74;q5=0.74;q6=0.74;  

Tya=0.1; td4=0.1;aw2=0.00;rw2=0.00;fp2=1000;


xxxxx0=0.01;yyyyy0=0.01;zzzzz0=0.01;wwwww0=0.01;ggggg0=0.01;mmmmm0=0.01;   %initial value

xxxxx(N+1)=[0];yyyyy(N+1)=[0];zzzzz(N+1)=[0];wwwww(N+1)=[0];ggggg(N+1)=[0];mmmmm(N+1)=[0];  %efficiency need improve
xxxxx1(N+1)=[0];yyyyy1(N+1)=[0];zzzzz1(N+1)=[0];wwwww1(N+1)=[0];ggggg1(N+1)=[0];mmmmm1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
xxxxx1(1)=xxxxx0+h^q1*(  yyyyy0      +vvvv1*xxxxx0+vvvv2*yyyyy0+vvvv3*zzzzz0+vvvv4*wwwww0+vvvv5*ggggg0+vvvv6*mmmmm0   )/(gamma(q1)*q1);
yyyyy1(1)=yyyyy0+h^q2*(  zzzzz0      +vvvv7*xxxxx0+vvvv8*yyyyy0+vvvv9*zzzzz0+vvvv10*wwwww0+vvvv11*ggggg0+vvvv12*mmmmm0   )/(gamma(q2)*q2);
zzzzz1(1)=zzzzz0+h^q3*(  -24*xxxxx0-24*yyyyy0-3*zzzzz0+mmmmm0        +vvvv13*xxxxx0+vvvv14*yyyyy0+vvvv15*zzzzz0+vvvv16*wwwww0+vvvv17*ggggg0+vvvv18*mmmmm0    )/(gamma(q3)*q3);
wwwww1(1)=wwwww0+h^q4*(  314*ggggg0    +vvvv19*xxxxx0+vvvv20*yyyyy0+vvvv21*zzzzz0+vvvv22*wwwww0+vvvv23*ggggg0+vvvv24*mmmmm0   )/(gamma(q4)*q4);
ggggg1(1)=ggggg0+h^q5*(  6.4*xxxxx0+0.8*zzzzz0-3/23*sin(wwwww0)+320/30429*sin(2*wwwww0)-2/9*ggggg0-7/45*mmmmm0    +vvvv25*xxxxx0+vvvv26*yyyyy0+vvvv27*zzzzz0+vvvv28*wwwww0+vvvv29*ggggg0+vvvv30*mmmmm0   )/(gamma(q5)*q5);
mmmmm1(1)=mmmmm0+h^q6*(  -1/Tya*mmmmm0   +vvvv31*xxxxx0+vvvv32*yyyyy0+vvvv33*zzzzz0+vvvv34*wwwww0+vvvv35*ggggg0+vvvv36*mmmmm0   )/(gamma(q6)*q6);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
xxxxx(1)=xxxxx0+h^q1*((  yyyyy0      +vvvv1*xxxxx0+vvvv2*yyyyy0+vvvv3*zzzzz0+vvvv4*wwwww0+vvvv5*ggggg0+vvvv6*mmmmm0  )+q1*(  yyyyy0      +vvvv1*xxxxx0+vvvv2*yyyyy0+vvvv3*zzzzz0+vvvv4*wwwww0+vvvv5*ggggg0+vvvv6*mmmmm0    ))/gamma(q1+2);
yyyyy(1)=yyyyy0+h^q2*((    zzzzz0      +vvvv7*xxxxx0+vvvv8*yyyyy0+vvvv9*zzzzz0+vvvv10*wwwww0+vvvv11*ggggg0+vvvv12*mmmmm0   )+q2*(  zzzzz0      +vvvv7*xxxxx0+vvvv8*yyyyy0+vvvv9*zzzzz0+vvvv10*wwwww0+vvvv11*ggggg0+vvvv12*mmmmm0     ))/gamma(q2+2);
zzzzz(1)=zzzzz0+h^q3*((    -24*xxxxx0-24*yyyyy0-3*zzzzz0+mmmmm0        +vvvv13*xxxxx0+vvvv14*yyyyy0+vvvv15*zzzzz0+vvvv16*wwwww0+vvvv17*ggggg0+vvvv18*mmmmm0   )+q3*(  -24*xxxxx0-24*yyyyy0-3*zzzzz0+mmmmm0        +vvvv13*xxxxx0+vvvv14*yyyyy0+vvvv15*zzzzz0+vvvv16*wwwww0+vvvv17*ggggg0+vvvv18*mmmmm0   ))/gamma(q3+2);
wwwww(1)=wwwww0+h^q4*((   314*ggggg0    +vvvv19*xxxxx0+vvvv20*yyyyy0+vvvv21*zzzzz0+vvvv22*wwwww0+vvvv23*ggggg0+vvvv24*mmmmm0   )+q4*(  314*ggggg0    +vvvv19*xxxxx0+vvvv20*yyyyy0+vvvv21*zzzzz0+vvvv22*wwwww0+vvvv23*ggggg0+vvvv24*mmmmm0   ))/gamma(q4+2);
ggggg(1)=ggggg0+h^q5*((   6.4*xxxxx0+0.8*zzzzz0-3/23*sin(wwwww0)+320/30429*sin(2*wwwww0)-2/9*ggggg0-7/45*mmmmm0  +vvvv25*xxxxx0+vvvv26*yyyyy0+vvvv27*zzzzz0+vvvv28*wwwww0+vvvv29*ggggg0+vvvv30*mmmmm0     )+q5*(  6.4*xxxxx0+0.8*zzzzz0-3/23*sin(wwwww0)+320/30429*sin(2*wwwww0)-2/9*ggggg0-7/45*mmmmm0  +vvvv25*xxxxx0+vvvv26*yyyyy0+vvvv27*zzzzz0+vvvv28*wwwww0+vvvv29*ggggg0+vvvv30*mmmmm0   ))/gamma(q5+2);
mmmmm(1)=mmmmm0+h^q6*((   -1/Tya*mmmmm0   +vvvv31*xxxxx0+vvvv32*yyyyy0+vvvv33*zzzzz0+vvvv34*wwwww0+vvvv35*ggggg0+vvvv36*mmmmm0  )+q6*(  -10*mmmmm0   +vvvv31*xxxxx0+vvvv32*yyyyy0+vvvv33*zzzzz0+vvvv34*wwwww0+vvvv35*ggggg0+vvvv36*mmmmm0   ))/gamma(q6+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td4/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(     yyyyy0      +vvvv1*xxxxx0+vvvv2*yyyyy0+vvvv3*zzzzz0+vvvv4*wwwww0+vvvv5*ggggg0+vvvv6*mmmmm0     );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     zzzzz0      +vvvv7*xxxxx0+vvvv8*yyyyy0+vvvv9*zzzzz0+vvvv10*wwwww0+vvvv11*ggggg0+vvvv12*mmmmm0    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -24*xxxxx0-24*yyyyy0-3*zzzzz0+mmmmm0        +vvvv13*xxxxx0+vvvv14*yyyyy0+vvvv15*zzzzz0+vvvv16*wwwww0+vvvv17*ggggg0+vvvv18*mmmmm0  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*ggggg0    +vvvv19*xxxxx0+vvvv20*yyyyy0+vvvv21*zzzzz0+vvvv22*wwwww0+vvvv23*ggggg0+vvvv24*mmmmm0  );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(    6.4*xxxxx0+0.8*zzzzz0-3/23*sin(wwwww0)+320/30429*sin(2*wwwww0)-2/9*ggggg0-7/45*mmmmm0  +vvvv25*xxxxx0+vvvv26*yyyyy0+vvvv27*zzzzz0+vvvv28*wwwww0+vvvv29*ggggg0+vvvv30*mmmmm0   );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(   -1/Tya*mmmmm0   +vvvv31*xxxxx0+vvvv32*yyyyy0+vvvv33*zzzzz0+vvvv34*wwwww0+vvvv35*ggggg0+vvvv36*mmmmm0    );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(    yyyyy0      +vvvv1*xxxxx0+vvvv2*yyyyy0+vvvv3*zzzzz0+vvvv4*wwwww0+vvvv5*ggggg0+vvvv6*mmmmm0  );
N2=((n+1)^q2-n^q2)*(     zzzzz0      +vvvv7*xxxxx0+vvvv8*yyyyy0+vvvv9*zzzzz0+vvvv10*wwwww0+vvvv11*ggggg0+vvvv12*mmmmm0     );
N3=((n+1)^q3-n^q3)*(    -24*xxxxx0-24*yyyyy0-3*zzzzz0+mmmmm0        +vvvv13*xxxxx0+vvvv14*yyyyy0+vvvv15*zzzzz0+vvvv16*wwwww0+vvvv17*ggggg0+vvvv18*mmmmm0 );
N4=((n+1)^q4-n^q4)*(      314*ggggg0    +vvvv19*xxxxx0+vvvv20*yyyyy0+vvvv21*zzzzz0+vvvv22*wwwww0+vvvv23*ggggg0+vvvv24*mmmmm0     );
N5=((n+1)^q5-n^q5)*(   6.4*xxxxx0+0.8*zzzzz0-3/23*sin(wwwww0)+320/30429*sin(2*wwwww0)-2/9*ggggg0-7/45*mmmmm0  +vvvv25*xxxxx0+vvvv26*yyyyy0+vvvv27*zzzzz0+vvvv28*wwwww0+vvvv29*ggggg0+vvvv30*mmmmm0   );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mmmmm0   +vvvv31*xxxxx0+vvvv32*yyyyy0+vvvv33*zzzzz0+vvvv34*wwwww0+vvvv35*ggggg0+vvvv36*mmmmm0       );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    yyyyy(j)      +vvvv1*xxxxx(j)+vvvv2*yyyyy(j)+vvvv3*zzzzz(j)+vvvv4*wwwww(j)+vvvv5*ggggg(j)+vvvv6*mmmmm(j)   );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     zzzzz(j)      +vvvv7*xxxxx(j)+vvvv8*yyyyy(j)+vvvv9*zzzzz(j)+vvvv10*wwwww(j)+vvvv11*ggggg(j)+vvvv12*mmmmm(j)     );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xxxxx(j)-24*yyyyy(j)-3*zzzzz(j)+mmmmm(j)        +vvvv13*xxxxx(j)+vvvv14*yyyyy(j)+vvvv15*zzzzz(j)+vvvv16*wwwww(j)+vvvv17*ggggg(j)+vvvv18*mmmmm(j)   );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+  314*ggggg(j)    +vvvv19*xxxxx(j)+vvvv20*yyyyy(j)+vvvv21*zzzzz(j)+vvvv22*wwwww(j)+vvvv23*ggggg(j)+vvvv24*mmmmm(j)   );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(  aw2*sin(fp2*j)+rw2*rand(1)+  6.4*xxxxx(j)+0.8*zzzzz(j)-3/23*sin(wwwww(j))+320/30429*sin(2*wwwww(j))-2/9*ggggg(j)-7/45*mmmmm(j)  +vvvv25*xxxxx(j)+vvvv26*yyyyy(j)+vvvv27*zzzzz(j)+vvvv28*wwwww(j)+vvvv29*ggggg(j)+vvvv30*mmmmm(j)      );
        M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(  aw2*sin(fp2*j)+rw2*rand(1)+-1/Tya*mmmmm(j)   +vvvv31*xxxxx(j)+vvvv32*yyyyy(j)+vvvv33*zzzzz(j)+vvvv34*wwwww(j)+vvvv35*ggggg(j)+vvvv36*mmmmm(j)    );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(      yyyyy(j)      +vvvv1*xxxxx(j)+vvvv2*yyyyy(j)+vvvv3*zzzzz(j)+vvvv4*wwwww(j)+vvvv5*ggggg(j)+vvvv6*mmmmm(j)    );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(    zzzzz(j)      +vvvv7*xxxxx(j)+vvvv8*yyyyy(j)+vvvv9*zzzzz(j)+vvvv10*wwwww(j)+vvvv11*ggggg(j)+vvvv12*mmmmm(j)  );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(    -24*xxxxx(j)-24*yyyyy(j)-3*zzzzz(j)+mmmmm(j)        +vvvv13*xxxxx(j)+vvvv14*yyyyy(j)+vvvv15*zzzzz(j)+vvvv16*wwwww(j)+vvvv17*ggggg(j)+vvvv18*mmmmm(j)    );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw2*sin(fp2*j)+rw2*rand(1)+  314*ggggg(j)    +vvvv19*xxxxx(j)+vvvv20*yyyyy(j)+vvvv21*zzzzz(j)+vvvv22*wwwww(j)+vvvv23*ggggg(j)+vvvv24*mmmmm(j)    );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw2*sin(fp2*j)+rw2*rand(1)+ 6.4*xxxxx(j)+0.8*zzzzz(j)-3/23*sin(wwwww(j))+320/30429*sin(2*wwwww(j))-2/9*ggggg(j)-7/45*mmmmm(j)  +vvvv25*xxxxx(j)+vvvv26*yyyyy(j)+vvvv27*zzzzz(j)+vvvv28*wwwww(j)+vvvv29*ggggg(j)+vvvv30*mmmmm(j)     );
        N6=N6+((n-j+1)^q6-(n-j)^q6)*(   aw2*sin(fp2*j)+rw2*rand(1)+ -1/Tya*mmmmm(j)   +vvvv31*xxxxx(j)+vvvv32*yyyyy(j)+vvvv33*zzzzz(j)+vvvv34*wwwww(j)+vvvv35*ggggg(j)+vvvv36*mmmmm(j)      );
end   
xxxxx1(n+1)=xxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
yyyyy1(n+1)=yyyyy0+h^q2*N2/(gamma(q2)*q2);
zzzzz1(n+1)=zzzzz0+h^q3*N3/(gamma(q3)*q3);
wwwww1(n+1)=wwwww0+h^q4*N4/(gamma(q4)*q4);
ggggg1(n+1)=ggggg0+h^q5*N5/(gamma(q5)*q5);
mmmmm1(n+1)=mmmmm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
xxxxx(n+2)=xxxxx0+h^q1*(   yyyyy(n+1)      +vvvv1*xxxxx(n+1)+vvvv2*yyyyy(n+1)+vvvv3*zzzzz(n+1)+vvvv4*wwwww(n+1)+vvvv5*ggggg(n+1)+vvvv6*mmmmm(n+1)   +M1)/gamma(q1+2);
yyyyy(n+2)=yyyyy0+h^q2*(    zzzzz(n+1)      +vvvv7*xxxxx(n+1)+vvvv8*yyyyy(n+1)+vvvv9*zzzzz(n+1)+vvvv10*wwwww(n+1)+vvvv11*ggggg(n+1)+vvvv12*mmmmm(n+1)        +M2)/gamma(q2+2);
zzzzz(n+2)=zzzzz0+h^q3*(     -24*xxxxx(n+1)-24*yyyyy(n+1)-3*zzzzz(n+1)+mmmmm(n+1)        +vvvv13*xxxxx(n+1)+vvvv14*yyyyy(n+1)+vvvv15*zzzzz(n+1)+vvvv16*wwwww(n+1)+vvvv17*ggggg(n+1)+vvvv18*mmmmm(n+1)    +M3)/gamma(q3+2);
wwwww(n+2)=wwwww0+h^q4*(     314*ggggg(n+1)    +vvvv19*xxxxx(n+1)+vvvv20*yyyyy(n+1)+vvvv21*zzzzz(n+1)+vvvv22*wwwww(n+1)+vvvv23*ggggg(n+1)+vvvv24*mmmmm(n+1)        +M4)/gamma(q4+2);
ggggg(n+2)=ggggg0+h^q5*(     6.4*xxxxx(n+1)+0.8*zzzzz(n+1)-3/23*sin(wwwww(n+1))+320/30429*sin(2*wwwww(n+1))-2/9*ggggg(n+1)-7/45*mmmmm(n+1)  +vvvv25*xxxxx(n+1)+vvvv26*yyyyy(n+1)+vvvv27*zzzzz(n+1)+vvvv28*wwwww(n+1)+vvvv29*ggggg(n+1)+vvvv30*mmmmm(n+1)    +M5)/gamma(q5+2);
mmmmm(n+2)=mmmmm0+h^q6*(     -1/Tya*mmmmm(n+1)   +vvvv31*xxxxx(n+1)+vvvv32*yyyyy(n+1)+vvvv33*zzzzz(n+1)+vvvv34*wwwww(n+1)+vvvv35*ggggg(n+1)+vvvv36*mmmmm(n+1)      +M6)/gamma(q6+2);

end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td4/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(    yyyyy0      +vvvv1*xxxxx0+vvvv2*yyyyy0+vvvv3*zzzzz0+vvvv4*wwwww0+vvvv5*ggggg0+vvvv6*mmmmm0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(      zzzzz0      +vvvv7*xxxxx0+vvvv8*yyyyy0+vvvv9*zzzzz0+vvvv10*wwwww0+vvvv11*ggggg0+vvvv12*mmmmm0     );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -24*xxxxx0-24*yyyyy0-3*zzzzz0+mmmmm0        +vvvv13*xxxxx0+vvvv14*yyyyy0+vvvv15*zzzzz0+vvvv16*wwwww0+vvvv17*ggggg0+vvvv18*mmmmm0   );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*ggggg0    +vvvv19*xxxxx0+vvvv20*yyyyy0+vvvv21*zzzzz0+vvvv22*wwwww0+vvvv23*ggggg0+vvvv24*mmmmm0   );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(     6.4*xxxxx0+0.8*zzzzz0-3/23*sin(wwwww0)+320/30429*sin(2*wwwww0)-2/9*ggggg0-7/45*mmmmm0  +vvvv25*xxxxx0+vvvv26*yyyyy0+vvvv27*zzzzz0+vvvv28*wwwww0+vvvv29*ggggg0+vvvv30*mmmmm0  );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(       -1/Tya*mmmmm0   +vvvv31*xxxxx0+vvvv32*yyyyy0+vvvv33*zzzzz0+vvvv34*wwwww0+vvvv35*ggggg0+vvvv36*mmmmm0    );

N1=((n+1)^q1-n^q1)*(     yyyyy0      +vvvv1*xxxxx0+vvvv2*yyyyy0+vvvv3*zzzzz0+vvvv4*wwwww0+vvvv5*ggggg0+vvvv6*mmmmm0       );
N2=((n+1)^q2-n^q2)*(      zzzzz0      +vvvv7*xxxxx0+vvvv8*yyyyy0+vvvv9*zzzzz0+vvvv10*wwwww0+vvvv11*ggggg0+vvvv12*mmmmm0      );
N3=((n+1)^q3-n^q3)*(     -24*xxxxx0-24*yyyyy0-3*zzzzz0+mmmmm0        +vvvv13*xxxxx0+vvvv14*yyyyy0+vvvv15*zzzzz0+vvvv16*wwwww0+vvvv17*ggggg0+vvvv18*mmmmm0  );
N4=((n+1)^q4-n^q4)*(    314*ggggg0    +vvvv19*xxxxx0+vvvv20*yyyyy0+vvvv21*zzzzz0+vvvv22*wwwww0+vvvv23*ggggg0+vvvv24*mmmmm0   );
N5=((n+1)^q5-n^q5)*(     6.4*xxxxx0+0.8*zzzzz0-3/23*sin(wwwww0)+320/30429*sin(2*wwwww0)-2/9*ggggg0-7/45*mmmmm0  +vvvv25*xxxxx0+vvvv26*yyyyy0+vvvv27*zzzzz0+vvvv28*wwwww0+vvvv29*ggggg0+vvvv30*mmmmm0  );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mmmmm0   +vvvv31*xxxxx0+vvvv32*yyyyy0+vvvv33*zzzzz0+vvvv34*wwwww0+vvvv35*ggggg0+vvvv36*mmmmm0  );
for j=1:td4/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace) 
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(      yyyyy(j)      +vvvv1*xxxxx(j)+vvvv2*yyyyy(j)+vvvv3*zzzzz(j)+vvvv4*wwwww(j)+vvvv5*ggggg(j)+vvvv6*mmmmm(j)      );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zzzzz(j)      +vvvv7*xxxxx(j)+vvvv8*yyyyy(j)+vvvv9*zzzzz(j)+vvvv10*wwwww(j)+vvvv11*ggggg(j)+vvvv12*mmmmm(j)  );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -24*xxxxx(j)-24*yyyyy(j)-3*zzzzz(j)+mmmmm(j)        +vvvv13*xxxxx(j)+vvvv14*yyyyy(j)+vvvv15*zzzzz(j)+vvvv16*wwwww(j)+vvvv17*ggggg(j)+vvvv18*mmmmm(j)  );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aw2*sin(fp2*j)+rw2*rand(1)+     314*ggggg(j)    +vvvv19*xxxxx(j)+vvvv20*yyyyy(j)+vvvv21*zzzzz(j)+vvvv22*wwwww(j)+vvvv23*ggggg(j)+vvvv24*mmmmm(j)     );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+  6.4*xxxxx(j)+0.8*zzzzz(j)-3/23*sin(wwwww(j))+320/30429*sin(2*wwwww(j))-2/9*ggggg(j)-7/45*mmmmm(j)  +vvvv25*xxxxx(j)+vvvv26*yyyyy(j)+vvvv27*zzzzz(j)+vvvv28*wwwww(j)+vvvv29*ggggg(j)+vvvv30*mmmmm(j)  );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+    -1/Tya*mmmmm(j)   +vvvv31*xxxxx(j)+vvvv32*yyyyy(j)+vvvv33*zzzzz(j)+vvvv34*wwwww(j)+vvvv35*ggggg(j)+vvvv36*mmmmm(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yyyyy(j)      +vvvv1*xxxxx(j)+vvvv2*yyyyy(j)+vvvv3*zzzzz(j)+vvvv4*wwwww(j)+vvvv5*ggggg(j)+vvvv6*mmmmm(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zzzzz(j)      +vvvv7*xxxxx(j)+vvvv8*yyyyy(j)+vvvv9*zzzzz(j)+vvvv10*wwwww(j)+vvvv11*ggggg(j)+vvvv12*mmmmm(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xxxxx(j)-24*yyyyy(j)-3*zzzzz(j)+mmmmm(j)        +vvvv13*xxxxx(j)+vvvv14*yyyyy(j)+vvvv15*zzzzz(j)+vvvv16*wwwww(j)+vvvv17*ggggg(j)+vvvv18*mmmmm(j)  );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw2*sin(fp2*j)+rw2*rand(1)+    314*ggggg(j)    +vvvv19*xxxxx(j)+vvvv20*yyyyy(j)+vvvv21*zzzzz(j)+vvvv22*wwwww(j)+vvvv23*ggggg(j)+vvvv24*mmmmm(j)       );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw2*sin(fp2*j)+rw2*rand(1)+    6.4*xxxxx(j)+0.8*zzzzz(j)-3/23*sin(wwwww(j))+320/30429*sin(2*wwwww(j))-2/9*ggggg(j)-7/45*mmmmm(j)  +vvvv25*xxxxx(j)+vvvv26*yyyyy(j)+vvvv27*zzzzz(j)+vvvv28*wwwww(j)+vvvv29*ggggg(j)+vvvv30*mmmmm(j) );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(   aw2*sin(fp2*j)+rw2*rand(1)+     -1/Tya*mmmmm(j)   +vvvv31*xxxxx(j)+vvvv32*yyyyy(j)+vvvv33*zzzzz(j)+vvvv34*wwwww(j)+vvvv35*ggggg(j)+vvvv36*mmmmm(j)       );
end   
for j=td4/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )       
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     yyyyy(j)      +vvvv1*xxxxx(j)+vvvv2*yyyyy(j)+vvvv3*zzzzz(j)+vvvv4*wwwww(j)+vvvv5*ggggg(j)+vvvv6*mmmmm(j)     );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zzzzz(j)      +vvvv7*xxxxx(j)+vvvv8*yyyyy(j)+vvvv9*zzzzz(j)+vvvv10*wwwww(j)+vvvv11*ggggg(j)+vvvv12*mmmmm(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xxxxx(j)-24*yyyyy(j)-3*zzzzz(j)+mmmmm(j-td4/h)        +vvvv13*xxxxx(j)+vvvv14*yyyyy(j)+vvvv15*zzzzz(j)+vvvv16*wwwww(j)+vvvv17*ggggg(j)+vvvv18*mmmmm(j) );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+   314*ggggg(j)    +vvvv19*xxxxx(j)+vvvv20*yyyyy(j)+vvvv21*zzzzz(j)+vvvv22*wwwww(j)+vvvv23*ggggg(j)+vvvv24*mmmmm(j)      );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(   aw2*sin(fp2*j)+rw2*rand(1)+    6.4*xxxxx(j)+0.8*zzzzz(j)-3/23*sin(wwwww(j))+320/30429*sin(2*wwwww(j))-2/9*ggggg(j)-7/45*mmmmm(j-td4/h)  +vvvv25*xxxxx(j)+vvvv26*yyyyy(j)+vvvv27*zzzzz(j)+vvvv28*wwwww(j)+vvvv29*ggggg(j)+vvvv30*mmmmm(j)   );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(    aw2*sin(fp2*j)+rw2*rand(1)+   -1/Tya*mmmmm(j-td4/h)   +vvvv31*xxxxx(j)+vvvv32*yyyyy(j)+vvvv33*zzzzz(j)+vvvv34*wwwww(j)+vvvv35*ggggg(j)+vvvv36*mmmmm(j)    );      
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yyyyy(j)      +vvvv1*xxxxx(j)+vvvv2*yyyyy(j)+vvvv3*zzzzz(j)+vvvv4*wwwww(j)+vvvv5*ggggg(j)+vvvv6*mmmmm(j)      );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zzzzz(j)      +vvvv7*xxxxx(j)+vvvv8*yyyyy(j)+vvvv9*zzzzz(j)+vvvv10*wwwww(j)+vvvv11*ggggg(j)+vvvv12*mmmmm(j)   );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xxxxx(j)-24*yyyyy(j)-3*zzzzz(j)+mmmmm(j-td4/h)        +vvvv13*xxxxx(j)+vvvv14*yyyyy(j)+vvvv15*zzzzz(j)+vvvv16*wwwww(j)+vvvv17*ggggg(j)+vvvv18*mmmmm(j)   );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw2*sin(fp2*j)+rw2*rand(1)+   314*ggggg(j)    +vvvv19*xxxxx(j)+vvvv20*yyyyy(j)+vvvv21*zzzzz(j)+vvvv22*wwwww(j)+vvvv23*ggggg(j)+vvvv24*mmmmm(j)   );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw2*sin(fp2*j)+rw2*rand(1)+   6.4*xxxxx(j)+0.8*zzzzz(j)-3/23*sin(wwwww(j))+320/30429*sin(2*wwwww(j))-2/9*ggggg(j)-7/45*mmmmm(j-td4/h)  +vvvv25*xxxxx(j)+vvvv26*yyyyy(j)+vvvv27*zzzzz(j)+vvvv28*wwwww(j)+vvvv29*ggggg(j)+vvvv30*mmmmm(j)  );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(    aw2*sin(fp2*j)+rw2*rand(1)+   -1/Tya*mmmmm(j-td4/h)   +vvvv31*xxxxx(j)+vvvv32*yyyyy(j)+vvvv33*zzzzz(j)+vvvv34*wwwww(j)+vvvv35*ggggg(j)+vvvv36*mmmmm(j)    );      
end
xxxxx1(n+1)=xxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
yyyyy1(n+1)=yyyyy0+h^q2*N2/(gamma(q2)*q2);
zzzzz1(n+1)=zzzzz0+h^q3*N3/(gamma(q2)*q3);
wwwww1(n+1)=wwwww0+h^q4*N4/(gamma(q2)*q4);
ggggg1(n+1)=ggggg0+h^q5*N5/(gamma(q5)*q5);
mmmmm1(n+1)=mmmmm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
xxxxx(n+1)=xxxxx0+h^q1*(     yyyyy(n+1)      +vvvv1*xxxxx(n+1)+vvvv2*yyyyy(n+1)+vvvv3*zzzzz(n+1)+vvvv4*wwwww(n+1)+vvvv5*ggggg(n+1)+vvvv6*mmmmm(n+1)   +M1)/gamma(q1+2);
yyyyy(n+1)=yyyyy0+h^q2*(     zzzzz(n+1)      +vvvv7*xxxxx(n+1)+vvvv8*yyyyy(n+1)+vvvv9*zzzzz(n+1)+vvvv10*wwwww(n+1)+vvvv11*ggggg(n+1)+vvvv12*mmmmm(n+1)     +M2)/gamma(q2+2);
zzzzz(n+1)=zzzzz0+h^q3*(     -24*xxxxx(n+1)-24*yyyyy(n+1)-3*zzzzz(n+1)+mmmmm(n+1-td4/h)        +vvvv13*xxxxx(n+1)+vvvv14*yyyyy(n+1)+vvvv15*zzzzz(n+1)+vvvv16*wwwww(n+1)+vvvv17*ggggg(n+1)+vvvv18*mmmmm(n+1)    +M3)/gamma(q3+2);
wwwww(n+1)=wwwww0+h^q4*(     314*ggggg(n+1)    +vvvv19*xxxxx(n+1)+vvvv20*yyyyy(n+1)+vvvv21*zzzzz(n+1)+vvvv22*wwwww(n+1)+vvvv23*ggggg(n+1)+vvvv24*mmmmm(n+1)      +M4)/gamma(q4+2);
ggggg(n+1)=ggggg0+h^q5*(     6.4*xxxxx(n+1)+0.8*zzzzz(n+1)-3/23*sin(wwwww(n+1))+320/30429*sin(2*wwwww(n+1))-2/9*ggggg(n+1)-7/45*mmmmm(n+1-td4/h)  +vvvv25*xxxxx(n+1)+vvvv26*yyyyy(n+1)+vvvv27*zzzzz(n+1)+vvvv28*wwwww(n+1)+vvvv29*ggggg(n+1)+vvvv30*mmmmm(n+1)    +M5)/gamma(q5+2);
mmmmm(n+1)=mmmmm0+h^q6*(     -1/Tya*mmmmm(n+1-td4/h)   +vvvv31*xxxxx(n+1)+vvvv32*yyyyy(n+1)+vvvv33*zzzzz(n+1)+vvvv34*wwwww(n+1)+vvvv35*ggggg(n+1)+vvvv36*mmmmm(n+1)      +M6)/gamma(q6+2);
end
%%
q1=1;q2=1;q3=1;q4=1;  
% q1=0.98;q2=0.98;q3=0.98;q4=0.98;  
aw=0.000386;rw=0.0001;fp=0.052;   %x6

h=0.01;N=999;           

px0=0.01;py0=0.01;pz0=0.01;pw0=0.01;   %initial value
td=0.1;   %0.08
kp=0.98;ki=0.8;kd=0.5;
a=314;Tab=19;D=2;Eq=1.35;xd=1.15;xq=1.474;
Tw=0.8;Ty=0.1;Vs=1;e1=0.5;e2=1;e=0.7;

px(N+1)=[0];py(N+1)=[0];pz(N+1)=[0];pw(N+1)=[0];  %efficiency need improve
px1(N+1)=[0];py1(N+1)=[0];pz1(N+1)=[0];pw1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model：
%  f( x0 , y0 , z0 delay term replace)
px1(1)=px0+h^q1*(           a*py0                                         )/(gamma(q1)*q1);
py1(1)=py0+h^q2*(           1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))         )/(gamma(q2)*q2);
pz1(1)=pz0+h^q3*(          1/(e1*Tw)*(-pz0+e2*pw0-(e*e2*Tw/Ty)*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0)))-pw0 )) )/(gamma(q3)*q3);
pw1(1)=pw0+h^q4*(          1/Ty*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))-pw0)) )/(gamma(q4)*q4);
% %    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
px(1)=px0+h^q1*((            a*py0              )+q1*(         a*py0              ))/gamma(q1+2);
py(1)=py0+h^q2*((           1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))     )+q2*(      1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))    ))/gamma(q2+2);
pz(1)=pz0+h^q3*((  1/(e1*Tw)*(-pz0+e2*pw0-(e*e2*Tw/Ty)*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0)))-pw0 )) )+q3*( 1/(e1*Tw)*(-pz0+e2*pw0-(e*e2*Tw/Ty)*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0)))-pw0 ))))/gamma(q3+2);
pw(1)=pw0+h^q4*((   1/Ty*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))-pw0))  )+q4*(   1/Ty*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))-pw0)) ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(         a*py0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(       1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pz0+e2*pw0-(e*e2*Tw/Ty)*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0)))-pw0 ))  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))-pw0)) );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(               a*py0     );
N2=((n+1)^q2-n^q2)*(     1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))     );
N3=((n+1)^q3-n^q3)*( 1/(e1*Tw)*(-pz0+e2*pw0-(e*e2*Tw/Ty)*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0)))-pw0 )));
N4=((n+1)^q4-n^q4)*(  1/Ty*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))-pw0)) );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(  aw*sin(fp*j)+rw*rand(1)+        a*py(j)        );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(  aw*sin(fp*j)+rw*rand(1)+   1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))   );
     M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pz(j)+e2*pw(j)-(e*e2*Tw/Ty)*(-kp*pw(j)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j))))-pw(j) )) );
     M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pw(j)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))-pw(j))) );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(   aw*sin(fp*j)+rw*rand(1)+            a*py(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(   aw*sin(fp*j)+rw*rand(1)+   1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))     );
    N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pz(j)+e2*pw(j)-(e*e2*Tw/Ty)*(-kp*pw(j)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j))))-pw(j) )) );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pw(j)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))-pw(j))));
end   
px1(n+1)=px0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
py1(n+1)=py0+h^q2*N2/(gamma(q2)*q2);
pz1(n+1)=pz0+h^q3*N3/(gamma(q3)*q3);
pw1(n+1)=pw0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

px(n+2)=px0+h^q1*(         a*py(n+1)    +M1)/gamma(q1+2);
py(n+2)=py0+h^q2*(       1/Tab*(pz(n+1)-D*py(n+1)-(Eq*Vs/xd)*sin(px(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(n+1)))   +M2)/gamma(q2+2);
pz(n+2)=pz0+h^q3*( 1/(e1*Tw)*(-pz(n+1)+e2*pw(n+1)-(e*e2*Tw/Ty)*(-kp*pw(n+1)-ki/a*px(n+1)-kd*(1/Tab*(pz(n+1)-D*py(n+1)-(Eq*Vs/xd)*sin(px(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(n+1))))-pw(n+1) )) +M3)/gamma(q3+2);
pw(n+2)=pw0+h^q4*( 1/Ty*(-kp*pw(n+1)-ki/a*px(n+1)-kd*(1/Tab*(pz(n+1)-D*py(n+1)-(Eq*Vs/xd)*sin(px(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(n+1)))-pw(n+1)))+M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(            a*py0           );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))   );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pz0+e2*pw0-(e*e2*Tw/Ty)*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0)))-pw0 )));
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))-pw0)));
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(           a*py0             );
N2=((n+1)^q2-n^q2)*(    1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))    );
N3=((n+1)^q3-n^q3)*(  1/(e1*Tw)*(-pz0+e2*pw0-(e*e2*Tw/Ty)*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0)))-pw0 )));
N4=((n+1)^q4-n^q4)*(   1/Ty*(-kp*pw0-ki/a*px0-kd*(1/Tab*(pz0-D*py0-(Eq*Vs/xd)*sin(px0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px0))-pw0)));
for j=1:td/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace)
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    aw*sin(fp*j)+rw*rand(1)+     a*py(j)               );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   aw*sin(fp*j)+rw*rand(1)+   1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))                    );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pz(j)+e2*pw(j)-(e*e2*Tw/Ty)*(-kp*pw(j)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j))))-pw(j) )) );
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pw(j)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))-pw(j))));
       %    f(  x(j) , y(j) , z(j)  delay term replace)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(   aw*sin(fp*j)+rw*rand(1)+     a*py(j)                  );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))      );
        N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pz(j)+e2*pw(j)-(e*e2*Tw/Ty)*(-kp*pw(j)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j))))-pw(j) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pw(j)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))-pw(j))));
     
end   
for j=td/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(   aw*sin(fp*j)+rw*rand(1)+      a*py(j)         );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))      );
        M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(aw*sin(fp*j)+rw*rand(1)+  1/(e1*Tw)*(-pz(j)+e2*pw(j-td/h)-(e*e2*Tw/Ty)*(-kp*pw(j-td/h)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j))))-pw(j-td/h) )));
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pw(j-td/h)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))-pw(j-td/h))) );
       %                   f(  x(j) , y(j) , z(j)  The delay variable is replaced by the value of the previous variable)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(  aw*sin(fp*j)+rw*rand(1)+      a*py(j)             );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))      );
       N3=N3+((n-j+1)^q3-(n-j)^q3)*( aw*sin(fp*j)+rw*rand(1)+1/(e1*Tw)*(-pz(j)+e2*pw(j-td/h)-(e*e2*Tw/Ty)*(-kp*pw(j-td/h)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j))))-pw(j-td/h) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pw(j-td/h)-ki/a*px(j)-kd*(1/Tab*(pz(j)-D*py(j)-(Eq*Vs/xd)*sin(px(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(j)))-pw(j-td/h))));
end

px1(n+1)=px0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
py1(n+1)=py0+h^q2*N2/(gamma(q2)*q2);
pz1(n+1)=pz0+h^q3*N3/(gamma(q3)*q3);
pw1(n+1)=pw0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
px(n+2)=px0+h^q1*(        a*py(n+1)            +M1)/gamma(q1+2);
py(n+2)=py0+h^q2*(     1/Tab*(pz(n+1)-D*py(n+1)-(Eq*Vs/xd)*sin(px(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(n+1)))   +M2)/gamma(q2+2);
pz(n+2)=pz0+h^q3*(  1/(e1*Tw)*(-pz(n+1)+e2*pw(n+1-td/h)-(e*e2*Tw/Ty)*(-kp*pw(n+1-td/h)-ki/a*px(n+1)-kd*(1/Tab*(pz(n+1)-D*py(n+1)-(Eq*Vs/xd)*sin(px(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(n+1))))-pw(n+1-td/h) ))+M3)/gamma(q3+2);
pw(n+2)=pw0+h^q4*(  1/Ty*(-kp*pw(n+1-td/h)-ki/a*px(n+1)-kd*(1/Tab*(pz(n+1)-D*py(n+1)-(Eq*Vs/xd)*sin(px(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*px(n+1)))-pw(n+1-td/h)))+M4)/gamma(q4+2);
end


q1=0.98;q2=0.98;q3=0.98;q4=0.98;  
% q1=0.95;q2=0.95;q3=0.95;q4=0.95; 
aw=0.000386;rw=0.0001;fp=0.052;   %x6previous variable

pxx0=0.01;pyy0=0.01;pzz0=0.01;pww0=0.01;   %initial value
td=0.1;   %0.08
kp=0.98;ki=0.8;kd=0.5;
a=314;Tab=19;D=2;Eq=1.35;xd=1.15;xq=1.474;
Tw=0.8;Ty=0.1;Vs=1;e1=0.5;e2=1;e=0.7;

pxx(N+1)=[0];pyy(N+1)=[0];pzz(N+1)=[0];pww(N+1)=[0];  %efficiency need improve
pxx1(N+1)=[0];pyy1(N+1)=[0];pzz1(N+1)=[0];pww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model：
%  f( x0 , y0 , z0 delay term replace)
pxx1(1)=pxx0+h^q1*(           a*pyy0                                         )/(gamma(q1)*q1);
pyy1(1)=pyy0+h^q2*(           1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))         )/(gamma(q2)*q2);
pzz1(1)=pzz0+h^q3*(          1/(e1*Tw)*(-pzz0+e2*pww0-(e*e2*Tw/Ty)*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0)))-pww0 )) )/(gamma(q3)*q3);
pww1(1)=pww0+h^q4*(          1/Ty*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))-pww0)) )/(gamma(q4)*q4);
% %    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
pxx(1)=pxx0+h^q1*((            a*pyy0              )+q1*(         a*pyy0              ))/gamma(q1+2);
pyy(1)=pyy0+h^q2*((           1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))     )+q2*(      1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))    ))/gamma(q2+2);
pzz(1)=pzz0+h^q3*((  1/(e1*Tw)*(-pzz0+e2*pww0-(e*e2*Tw/Ty)*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0)))-pww0 )) )+q3*( 1/(e1*Tw)*(-pzz0+e2*pww0-(e*e2*Tw/Ty)*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0)))-pww0 ))))/gamma(q3+2);
pww(1)=pww0+h^q4*((   1/Ty*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))-pww0))  )+q4*(   1/Ty*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))-pww0)) ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(         a*pyy0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(       1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pzz0+e2*pww0-(e*e2*Tw/Ty)*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0)))-pww0 ))  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))-pww0)) );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(               a*pyy0     );
N2=((n+1)^q2-n^q2)*(     1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))     );
N3=((n+1)^q3-n^q3)*( 1/(e1*Tw)*(-pzz0+e2*pww0-(e*e2*Tw/Ty)*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0)))-pww0 )));
N4=((n+1)^q4-n^q4)*(  1/Ty*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))-pww0)) );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(  aw*sin(fp*j)+rw*rand(1)+        a*pyy(j)        );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*( aw*sin(fp*j)+rw*rand(1)+    1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))   );
     M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzz(j)+e2*pww(j)-(e*e2*Tw/Ty)*(-kp*pww(j)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j))))-pww(j) )) );
     M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pww(j)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))-pww(j))) );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(   aw*sin(fp*j)+rw*rand(1)+            a*pyy(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(   aw*sin(fp*j)+rw*rand(1)+   1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))     );
    N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzz(j)+e2*pww(j)-(e*e2*Tw/Ty)*(-kp*pww(j)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j))))-pww(j) )) );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pww(j)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))-pww(j))));
end   
pxx1(n+1)=pxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
pyy1(n+1)=pyy0+h^q2*N2/(gamma(q2)*q2);
pzz1(n+1)=pzz0+h^q3*N3/(gamma(q3)*q3);
pww1(n+1)=pww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

pxx(n+2)=pxx0+h^q1*(         a*pyy(n+1)    +M1)/gamma(q1+2);
pyy(n+2)=pyy0+h^q2*(       1/Tab*(pzz(n+1)-D*pyy(n+1)-(Eq*Vs/xd)*sin(pxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(n+1)))   +M2)/gamma(q2+2);
pzz(n+2)=pzz0+h^q3*( 1/(e1*Tw)*(-pzz(n+1)+e2*pww(n+1)-(e*e2*Tw/Ty)*(-kp*pww(n+1)-ki/a*pxx(n+1)-kd*(1/Tab*(pzz(n+1)-D*pyy(n+1)-(Eq*Vs/xd)*sin(pxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(n+1))))-pww(n+1) )) +M3)/gamma(q3+2);
pww(n+2)=pww0+h^q4*( 1/Ty*(-kp*pww(n+1)-ki/a*pxx(n+1)-kd*(1/Tab*(pzz(n+1)-D*pyy(n+1)-(Eq*Vs/xd)*sin(pxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(n+1)))-pww(n+1)))+M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(            a*pyy0           );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))   );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pzz0+e2*pww0-(e*e2*Tw/Ty)*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0)))-pww0 )));
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))-pww0)));
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(           a*pyy0             );
N2=((n+1)^q2-n^q2)*(    1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))    );
N3=((n+1)^q3-n^q3)*(  1/(e1*Tw)*(-pzz0+e2*pww0-(e*e2*Tw/Ty)*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0)))-pww0 )));
N4=((n+1)^q4-n^q4)*(   1/Ty*(-kp*pww0-ki/a*pxx0-kd*(1/Tab*(pzz0-D*pyy0-(Eq*Vs/xd)*sin(pxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx0))-pww0)));
for j=1:td/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace)
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(  aw*sin(fp*j)+rw*rand(1)+       a*pyy(j)               );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(  aw*sin(fp*j)+rw*rand(1)+    1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))                    );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzz(j)+e2*pww(j)-(e*e2*Tw/Ty)*(-kp*pww(j)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j))))-pww(j) )) );
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pww(j)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))-pww(j))));
       %    f(  x(j) , y(j) , z(j)  delay term replace)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(   aw*sin(fp*j)+rw*rand(1)+     a*pyy(j)                  );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))      );
        N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzz(j)+e2*pww(j)-(e*e2*Tw/Ty)*(-kp*pww(j)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j))))-pww(j) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pww(j)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))-pww(j))));
     
end   
for j=td/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(  aw*sin(fp*j)+rw*rand(1)+       a*pyy(j)         );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(  aw*sin(fp*j)+rw*rand(1)+   1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))      );
        M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzz(j)+e2*pww(j-td/h)-(e*e2*Tw/Ty)*(-kp*pww(j-td/h)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j))))-pww(j-td/h) )));
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pww(j-td/h)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))-pww(j-td/h))) );
       %                   f(  x(j) , y(j) , z(j)  The delay variable is replaced by the value of the previous variable)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(  aw*sin(fp*j)+rw*rand(1)+      a*pyy(j)             );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))      );
       N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzz(j)+e2*pww(j-td/h)-(e*e2*Tw/Ty)*(-kp*pww(j-td/h)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j))))-pww(j-td/h) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*(aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pww(j-td/h)-ki/a*pxx(j)-kd*(1/Tab*(pzz(j)-D*pyy(j)-(Eq*Vs/xd)*sin(pxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(j)))-pww(j-td/h))));
end

pxx1(n+1)=pxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
pyy1(n+1)=pyy0+h^q2*N2/(gamma(q2)*q2);
pzz1(n+1)=pzz0+h^q3*N3/(gamma(q3)*q3);
pww1(n+1)=pww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
pxx(n+2)=pxx0+h^q1*(        a*pyy(n+1)            +M1)/gamma(q1+2);
pyy(n+2)=pyy0+h^q2*(     1/Tab*(pzz(n+1)-D*pyy(n+1)-(Eq*Vs/xd)*sin(pxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(n+1)))   +M2)/gamma(q2+2);
pzz(n+2)=pzz0+h^q3*(  1/(e1*Tw)*(-pzz(n+1)+e2*pww(n+1-td/h)-(e*e2*Tw/Ty)*(-kp*pww(n+1-td/h)-ki/a*pxx(n+1)-kd*(1/Tab*(pzz(n+1)-D*pyy(n+1)-(Eq*Vs/xd)*sin(pxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(n+1))))-pww(n+1-td/h) ))+M3)/gamma(q3+2);
pww(n+2)=pww0+h^q4*(  1/Ty*(-kp*pww(n+1-td/h)-ki/a*pxx(n+1)-kd*(1/Tab*(pzz(n+1)-D*pyy(n+1)-(Eq*Vs/xd)*sin(pxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxx(n+1)))-pww(n+1-td/h)))+M4)/gamma(q4+2);
end

q1=0.95;q2=0.95;q3=0.95;q4=0.95;  
% q1=0.9;q2=0.9;q3=0.9;q4=0.9;  
aw=0.0003;rw=0.0001;fp=0.052;   %x6   0.7

pxxx0=0.01;pyyy0=0.01;pzzz0=0.01;pwww0=0.01;   %initial value
td=0.15;   %0.08
kp=0.98;ki=0.8;kd=0.5;
a=314;Tab=19;D=2;Eq=1.35;xd=1.15;xq=1.474;
Tw=0.8;Ty=0.1;Vs=1;e1=0.5;e2=1;e=0.7;

pxxx(N+1)=[0];pyyy(N+1)=[0];pzzz(N+1)=[0];pwww(N+1)=[0];  %efficiency need improve
pxxx1(N+1)=[0];pyyy1(N+1)=[0];pzzz1(N+1)=[0];pwww1(N+1)=[0];


%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model：
%  f( x0 , y0 , z0 delay term replace)
pxxx1(1)=pxxx0+h^q1*(           a*pyyy0                                         )/(gamma(q1)*q1);
pyyy1(1)=pyyy0+h^q2*(           1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))         )/(gamma(q2)*q2);
pzzz1(1)=pzzz0+h^q3*(          1/(e1*Tw)*(-pzzz0+e2*pwww0-(e*e2*Tw/Ty)*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0)))-pwww0 )) )/(gamma(q3)*q3);
pwww1(1)=pwww0+h^q4*(          1/Ty*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))-pwww0)) )/(gamma(q4)*q4);
% %    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
pxxx(1)=pxxx0+h^q1*((            a*pyyy0              )+q1*(         a*pyyy0              ))/gamma(q1+2);
pyyy(1)=pyyy0+h^q2*((           1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))     )+q2*(      1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))    ))/gamma(q2+2);
pzzz(1)=pzzz0+h^q3*((  1/(e1*Tw)*(-pzzz0+e2*pwww0-(e*e2*Tw/Ty)*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0)))-pwww0 )) )+q3*( 1/(e1*Tw)*(-pzzz0+e2*pwww0-(e*e2*Tw/Ty)*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0)))-pwww0 ))))/gamma(q3+2);
pwww(1)=pwww0+h^q4*((   1/Ty*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))-pwww0))  )+q4*(   1/Ty*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))-pwww0)) ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(         a*pyyy0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(       1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pzzz0+e2*pwww0-(e*e2*Tw/Ty)*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0)))-pwww0 ))  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))-pwww0)) );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(               a*pyyy0     );
N2=((n+1)^q2-n^q2)*(     1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))     );
N3=((n+1)^q3-n^q3)*( 1/(e1*Tw)*(-pzzz0+e2*pwww0-(e*e2*Tw/Ty)*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0)))-pwww0 )));
N4=((n+1)^q4-n^q4)*(  1/Ty*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))-pwww0)) );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(   aw*sin(fp*j)+rw*rand(1)+       a*pyyy(j)        );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*( aw*sin(fp*j)+rw*rand(1)+    1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))   );
     M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(aw*sin(fp*j)+rw*rand(1)+  1/(e1*Tw)*(-pzzz(j)+e2*pwww(j)-(e*e2*Tw/Ty)*(-kp*pwww(j)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j))))-pwww(j) )) );
     M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pwww(j)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))-pwww(j))) );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(  aw*sin(fp*j)+rw*rand(1)+             a*pyyy(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+    1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))     );
    N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzz(j)+e2*pwww(j)-(e*e2*Tw/Ty)*(-kp*pwww(j)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j))))-pwww(j) )) );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwww(j)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))-pwww(j))));
end   
pxxx1(n+1)=pxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
pyyy1(n+1)=pyyy0+h^q2*N2/(gamma(q2)*q2);
pzzz1(n+1)=pzzz0+h^q3*N3/(gamma(q3)*q3);
pwww1(n+1)=pwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

pxxx(n+2)=pxxx0+h^q1*(         a*pyyy(n+1)    +M1)/gamma(q1+2);
pyyy(n+2)=pyyy0+h^q2*(       1/Tab*(pzzz(n+1)-D*pyyy(n+1)-(Eq*Vs/xd)*sin(pxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(n+1)))   +M2)/gamma(q2+2);
pzzz(n+2)=pzzz0+h^q3*( 1/(e1*Tw)*(-pzzz(n+1)+e2*pwww(n+1)-(e*e2*Tw/Ty)*(-kp*pwww(n+1)-ki/a*pxxx(n+1)-kd*(1/Tab*(pzzz(n+1)-D*pyyy(n+1)-(Eq*Vs/xd)*sin(pxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(n+1))))-pwww(n+1) )) +M3)/gamma(q3+2);
pwww(n+2)=pwww0+h^q4*( 1/Ty*(-kp*pwww(n+1)-ki/a*pxxx(n+1)-kd*(1/Tab*(pzzz(n+1)-D*pyyy(n+1)-(Eq*Vs/xd)*sin(pxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(n+1)))-pwww(n+1)))+M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(            a*pyyy0           );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))   );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pzzz0+e2*pwww0-(e*e2*Tw/Ty)*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0)))-pwww0 )));
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))-pwww0)));
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(           a*pyyy0             );
N2=((n+1)^q2-n^q2)*(    1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))    );
N3=((n+1)^q3-n^q3)*(  1/(e1*Tw)*(-pzzz0+e2*pwww0-(e*e2*Tw/Ty)*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0)))-pwww0 )));
N4=((n+1)^q4-n^q4)*(   1/Ty*(-kp*pwww0-ki/a*pxxx0-kd*(1/Tab*(pzzz0-D*pyyy0-(Eq*Vs/xd)*sin(pxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx0))-pwww0)));
for j=1:td/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace)
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*( aw*sin(fp*j)+rw*rand(1)+        a*pyyy(j)               );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(  aw*sin(fp*j)+rw*rand(1)+    1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))                    );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzz(j)+e2*pwww(j)-(e*e2*Tw/Ty)*(-kp*pwww(j)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j))))-pwww(j) )) );
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwww(j)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))-pwww(j))));
       %    f(  x(j) , y(j) , z(j)  delay term replace)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(  aw*sin(fp*j)+rw*rand(1)+      a*pyyy(j)                  );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(   aw*sin(fp*j)+rw*rand(1)+ 1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))      );
        N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzz(j)+e2*pwww(j)-(e*e2*Tw/Ty)*(-kp*pwww(j)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j))))-pwww(j) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwww(j)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))-pwww(j))));
     
end   
for j=td/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(   aw*sin(fp*j)+rw*rand(1)+      a*pyyy(j)         );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))      );
        M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(aw*sin(fp*j)+rw*rand(1)+  1/(e1*Tw)*(-pzzz(j)+e2*pwww(j-td/h)-(e*e2*Tw/Ty)*(-kp*pwww(j-td/h)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j))))-pwww(j-td/h) )));
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwww(j-td/h)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))-pwww(j-td/h))) );
       %                   f(  x(j) , y(j) , z(j)  The delay variable is replaced by the value of the previous variable)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(  aw*sin(fp*j)+rw*rand(1)+      a*pyyy(j)             );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))      );
       N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzz(j)+e2*pwww(j-td/h)-(e*e2*Tw/Ty)*(-kp*pwww(j-td/h)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j))))-pwww(j-td/h) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*(aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pwww(j-td/h)-ki/a*pxxx(j)-kd*(1/Tab*(pzzz(j)-D*pyyy(j)-(Eq*Vs/xd)*sin(pxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(j)))-pwww(j-td/h))));
end

pxxx1(n+1)=pxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
pyyy1(n+1)=pyyy0+h^q2*N2/(gamma(q2)*q2);
pzzz1(n+1)=pzzz0+h^q3*N3/(gamma(q3)*q3);
pwww1(n+1)=pwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
pxxx(n+2)=pxxx0+h^q1*(        a*pyyy(n+1)            +M1)/gamma(q1+2);
pyyy(n+2)=pyyy0+h^q2*(     1/Tab*(pzzz(n+1)-D*pyyy(n+1)-(Eq*Vs/xd)*sin(pxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(n+1)))   +M2)/gamma(q2+2);
pzzz(n+2)=pzzz0+h^q3*(  1/(e1*Tw)*(-pzzz(n+1)+e2*pwww(n+1-td/h)-(e*e2*Tw/Ty)*(-kp*pwww(n+1-td/h)-ki/a*pxxx(n+1)-kd*(1/Tab*(pzzz(n+1)-D*pyyy(n+1)-(Eq*Vs/xd)*sin(pxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(n+1))))-pwww(n+1-td/h) ))+M3)/gamma(q3+2);
pwww(n+2)=pwww0+h^q4*(  1/Ty*(-kp*pwww(n+1-td/h)-ki/a*pxxx(n+1)-kd*(1/Tab*(pzzz(n+1)-D*pyyy(n+1)-(Eq*Vs/xd)*sin(pxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxx(n+1)))-pwww(n+1-td/h)))+M4)/gamma(q4+2);
end

q1=0.92;q2=0.92;q3=0.92;q4=0.92;  
% q1=0.89;q2=0.89;q3=0.89;q4=0.89;  
aw=0.00038;rw=0.0001;fp=0.052;   %x6previous variable    0.6iable

pxxxx0=0.01;pyyyy0=0.01;pzzzz0=0.01;pwwww0=0.01;   %initial value
td=0.13;   %0.08
kp=0.98;ki=0.8;kd=0.5;
a=314;Tab=19;D=2;Eq=1.35;xd=1.15;xq=1.474;
Tw=0.8;Ty=0.1;Vs=1;e1=0.5;e2=1;e=0.7;

pxxxx(N+1)=[0];pyyyy(N+1)=[0];pzzzz(N+1)=[0];pwwww(N+1)=[0];  %efficiency need improve
pxxxx1(N+1)=[0];pyyyy1(N+1)=[0];pzzzz1(N+1)=[0];pwwww1(N+1)=[0];


%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model：
%  f( x0 , y0 , z0 delay term replace)
pxxxx1(1)=pxxxx0+h^q1*(           a*pyyyy0                                         )/(gamma(q1)*q1);
pyyyy1(1)=pyyyy0+h^q2*(           1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))         )/(gamma(q2)*q2);
pzzzz1(1)=pzzzz0+h^q3*(          1/(e1*Tw)*(-pzzzz0+e2*pwwww0-(e*e2*Tw/Ty)*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0)))-pwwww0 )) )/(gamma(q3)*q3);
pwwww1(1)=pwwww0+h^q4*(          1/Ty*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))-pwwww0)) )/(gamma(q4)*q4);
% %    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
pxxxx(1)=pxxxx0+h^q1*((            a*pyyyy0              )+q1*(         a*pyyyy0              ))/gamma(q1+2);
pyyyy(1)=pyyyy0+h^q2*((           1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))     )+q2*(      1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))    ))/gamma(q2+2);
pzzzz(1)=pzzzz0+h^q3*((  1/(e1*Tw)*(-pzzzz0+e2*pwwww0-(e*e2*Tw/Ty)*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0)))-pwwww0 )) )+q3*( 1/(e1*Tw)*(-pzzzz0+e2*pwwww0-(e*e2*Tw/Ty)*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0)))-pwwww0 ))))/gamma(q3+2);
pwwww(1)=pwwww0+h^q4*((   1/Ty*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))-pwwww0))  )+q4*(   1/Ty*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))-pwwww0)) ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(         a*pyyyy0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(       1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pzzzz0+e2*pwwww0-(e*e2*Tw/Ty)*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0)))-pwwww0 ))  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))-pwwww0)) );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(               a*pyyyy0     );
N2=((n+1)^q2-n^q2)*(     1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))     );
N3=((n+1)^q3-n^q3)*( 1/(e1*Tw)*(-pzzzz0+e2*pwwww0-(e*e2*Tw/Ty)*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0)))-pwwww0 )));
N4=((n+1)^q4-n^q4)*(  1/Ty*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))-pwwww0)) );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(   aw*sin(fp*j)+rw*rand(1)+       a*pyyyy(j)        );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*( aw*sin(fp*j)+rw*rand(1)+    1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))   );
     M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzz(j)+e2*pwwww(j)-(e*e2*Tw/Ty)*(-kp*pwwww(j)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j))))-pwwww(j) )) );
     M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pwwww(j)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))-pwwww(j))) );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(  aw*sin(fp*j)+rw*rand(1)+             a*pyyyy(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+    1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))     );
    N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzz(j)+e2*pwwww(j)-(e*e2*Tw/Ty)*(-kp*pwwww(j)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j))))-pwwww(j) )) );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwwww(j)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))-pwwww(j))));
end   
pxxxx1(n+1)=pxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
pyyyy1(n+1)=pyyyy0+h^q2*N2/(gamma(q2)*q2);
pzzzz1(n+1)=pzzzz0+h^q3*N3/(gamma(q3)*q3);
pwwww1(n+1)=pwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

pxxxx(n+2)=pxxxx0+h^q1*(         a*pyyyy(n+1)    +M1)/gamma(q1+2);
pyyyy(n+2)=pyyyy0+h^q2*(       1/Tab*(pzzzz(n+1)-D*pyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(n+1)))   +M2)/gamma(q2+2);
pzzzz(n+2)=pzzzz0+h^q3*( 1/(e1*Tw)*(-pzzzz(n+1)+e2*pwwww(n+1)-(e*e2*Tw/Ty)*(-kp*pwwww(n+1)-ki/a*pxxxx(n+1)-kd*(1/Tab*(pzzzz(n+1)-D*pyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(n+1))))-pwwww(n+1) )) +M3)/gamma(q3+2);
pwwww(n+2)=pwwww0+h^q4*( 1/Ty*(-kp*pwwww(n+1)-ki/a*pxxxx(n+1)-kd*(1/Tab*(pzzzz(n+1)-D*pyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(n+1)))-pwwww(n+1)))+M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(            a*pyyyy0           );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))   );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pzzzz0+e2*pwwww0-(e*e2*Tw/Ty)*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0)))-pwwww0 )));
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))-pwwww0)));
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(           a*pyyyy0             );
N2=((n+1)^q2-n^q2)*(    1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))    );
N3=((n+1)^q3-n^q3)*(  1/(e1*Tw)*(-pzzzz0+e2*pwwww0-(e*e2*Tw/Ty)*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0)))-pwwww0 )));
N4=((n+1)^q4-n^q4)*(   1/Ty*(-kp*pwwww0-ki/a*pxxxx0-kd*(1/Tab*(pzzzz0-D*pyyyy0-(Eq*Vs/xd)*sin(pxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx0))-pwwww0)));
for j=1:td/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace)
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(   aw*sin(fp*j)+rw*rand(1)+      a*pyyyy(j)               );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   aw*sin(fp*j)+rw*rand(1)+   1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))                    );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzz(j)+e2*pwwww(j)-(e*e2*Tw/Ty)*(-kp*pwwww(j)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j))))-pwwww(j) )) );
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwwww(j)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))-pwwww(j))));
       %    f(  x(j) , y(j) , z(j)  delay term replace)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(   aw*sin(fp*j)+rw*rand(1)+     a*pyyyy(j)                  );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))      );
        N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzz(j)+e2*pwwww(j)-(e*e2*Tw/Ty)*(-kp*pwwww(j)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j))))-pwwww(j) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwwww(j)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))-pwwww(j))));
     
end   
for j=td/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(  aw*sin(fp*j)+rw*rand(1)+       a*pyyyy(j)         );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(  aw*sin(fp*j)+rw*rand(1)+   1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))      );
        M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzz(j)+e2*pwwww(j-td/h)-(e*e2*Tw/Ty)*(-kp*pwwww(j-td/h)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j))))-pwwww(j-td/h) )));
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pwwww(j-td/h)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))-pwwww(j-td/h))) );
       %                   f(  x(j) , y(j) , z(j)  The delay variable is replaced by the value of the previous variable)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(  aw*sin(fp*j)+rw*rand(1)+      a*pyyyy(j)             );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))      );
       N3=N3+((n-j+1)^q3-(n-j)^q3)*( aw*sin(fp*j)+rw*rand(1)+1/(e1*Tw)*(-pzzzz(j)+e2*pwwww(j-td/h)-(e*e2*Tw/Ty)*(-kp*pwwww(j-td/h)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j))))-pwwww(j-td/h) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwwww(j-td/h)-ki/a*pxxxx(j)-kd*(1/Tab*(pzzzz(j)-D*pyyyy(j)-(Eq*Vs/xd)*sin(pxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(j)))-pwwww(j-td/h))));
end

pxxxx1(n+1)=pxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
pyyyy1(n+1)=pyyyy0+h^q2*N2/(gamma(q2)*q2);
pzzzz1(n+1)=pzzzz0+h^q3*N3/(gamma(q3)*q3);
pwwww1(n+1)=pwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
pxxxx(n+2)=pxxxx0+h^q1*(        a*pyyyy(n+1)            +M1)/gamma(q1+2);
pyyyy(n+2)=pyyyy0+h^q2*(     1/Tab*(pzzzz(n+1)-D*pyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(n+1)))   +M2)/gamma(q2+2);
pzzzz(n+2)=pzzzz0+h^q3*(  1/(e1*Tw)*(-pzzzz(n+1)+e2*pwwww(n+1-td/h)-(e*e2*Tw/Ty)*(-kp*pwwww(n+1-td/h)-ki/a*pxxxx(n+1)-kd*(1/Tab*(pzzzz(n+1)-D*pyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(n+1))))-pwwww(n+1-td/h) ))+M3)/gamma(q3+2);
pwwww(n+2)=pwwww0+h^q4*(  1/Ty*(-kp*pwwww(n+1-td/h)-ki/a*pxxxx(n+1)-kd*(1/Tab*(pzzzz(n+1)-D*pyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxx(n+1)))-pwwww(n+1-td/h)))+M4)/gamma(q4+2);
end

% q1=0.89;q2=0.89;q3=0.89;q4=0.89;  
q1=0.88;q2=0.88;q3=0.88;q4=0.88;  
aw=0.000386;rw=0.0001;fp=0.052;   %x6previous variable

pxxxxx0=0.01;pyyyyy0=0.01;pzzzzz0=0.01;pwwwww0=0.01;   %initial value
td=0.13;   %0.08
kp=0.98;ki=0.8;kd=0.5;
a=314;Tab=19;D=2;Eq=1.35;xd=1.15;xq=1.474;
Tw=0.8;Ty=0.1;Vs=1;e1=0.5;e2=1;e=0.7;

pxxxxx(N+1)=[0];pyyyyy(N+1)=[0];pzzzzz(N+1)=[0];pwwwww(N+1)=[0];  %efficiency need improve
pxxxxx1(N+1)=[0];pyyyyy1(N+1)=[0];pzzzzz1(N+1)=[0];pwwwww1(N+1)=[0];


%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model：
%  f( x0 , y0 , z0 delay term replace)
pxxxxx1(1)=pxxxxx0+h^q1*(           a*pyyyyy0                                         )/(gamma(q1)*q1);
pyyyyy1(1)=pyyyyy0+h^q2*(           1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))         )/(gamma(q2)*q2);
pzzzzz1(1)=pzzzzz0+h^q3*(          1/(e1*Tw)*(-pzzzzz0+e2*pwwwww0-(e*e2*Tw/Ty)*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0)))-pwwwww0 )) )/(gamma(q3)*q3);
pwwwww1(1)=pwwwww0+h^q4*(          1/Ty*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))-pwwwww0)) )/(gamma(q4)*q4);
% %    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
pxxxxx(1)=pxxxxx0+h^q1*((            a*pyyyyy0              )+q1*(         a*pyyyyy0              ))/gamma(q1+2);
pyyyyy(1)=pyyyyy0+h^q2*((           1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))     )+q2*(      1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))    ))/gamma(q2+2);
pzzzzz(1)=pzzzzz0+h^q3*((  1/(e1*Tw)*(-pzzzzz0+e2*pwwwww0-(e*e2*Tw/Ty)*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0)))-pwwwww0 )) )+q3*( 1/(e1*Tw)*(-pzzzzz0+e2*pwwwww0-(e*e2*Tw/Ty)*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0)))-pwwwww0 ))))/gamma(q3+2);
pwwwww(1)=pwwwww0+h^q4*((   1/Ty*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))-pwwwww0))  )+q4*(   1/Ty*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))-pwwwww0)) ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value ++++++
for n=0:td/h-1               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(         a*pyyyyy0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(       1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pzzzzz0+e2*pwwwww0-(e*e2*Tw/Ty)*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0)))-pwwwww0 ))  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))-pwwwww0)) );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(               a*pyyyyy0     );
N2=((n+1)^q2-n^q2)*(     1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))     );
N3=((n+1)^q3-n^q3)*( 1/(e1*Tw)*(-pzzzzz0+e2*pwwwww0-(e*e2*Tw/Ty)*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0)))-pwwwww0 )));
N4=((n+1)^q4-n^q4)*(  1/Ty*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))-pwwwww0)) );
for j=1:n   %    f(  x(j) , y(j) , z(j) delay term replace )
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*( aw*sin(fp*j)+rw*rand(1)+         a*pyyyyy(j)        );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(  aw*sin(fp*j)+rw*rand(1)+   1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))   );
     M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzzz(j)+e2*pwwwww(j)-(e*e2*Tw/Ty)*(-kp*pwwwww(j)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j))))-pwwwww(j) )) );
     M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pwwwww(j)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))-pwwwww(j))) );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(  aw*sin(fp*j)+rw*rand(1)+             a*pyyyyy(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+    1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))     );
    N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzzz(j)+e2*pwwwww(j)-(e*e2*Tw/Ty)*(-kp*pwwwww(j)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j))))-pwwwww(j) )) );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwwwww(j)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))-pwwwww(j))));
end   
pxxxxx1(n+1)=pxxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
pyyyyy1(n+1)=pyyyyy0+h^q2*N2/(gamma(q2)*q2);
pzzzzz1(n+1)=pzzzzz0+h^q3*N3/(gamma(q3)*q3);
pwwwww1(n+1)=pwwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

pxxxxx(n+2)=pxxxxx0+h^q1*(         a*pyyyyy(n+1)    +M1)/gamma(q1+2);
pyyyyy(n+2)=pyyyyy0+h^q2*(       1/Tab*(pzzzzz(n+1)-D*pyyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(n+1)))   +M2)/gamma(q2+2);
pzzzzz(n+2)=pzzzzz0+h^q3*( 1/(e1*Tw)*(-pzzzzz(n+1)+e2*pwwwww(n+1)-(e*e2*Tw/Ty)*(-kp*pwwwww(n+1)-ki/a*pxxxxx(n+1)-kd*(1/Tab*(pzzzzz(n+1)-D*pyyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(n+1))))-pwwwww(n+1) )) +M3)/gamma(q3+2);
pwwwww(n+2)=pwwwww0+h^q4*( 1/Ty*(-kp*pwwwww(n+1)-ki/a*pxxxxx(n+1)-kd*(1/Tab*(pzzzzz(n+1)-D*pyyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(n+1)))-pwwwww(n+1)))+M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------
for n=td/h:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(            a*pyyyyy0           );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))   );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*( 1/(e1*Tw)*(-pzzzzz0+e2*pwwwww0-(e*e2*Tw/Ty)*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0)))-pwwwww0 )));
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   1/Ty*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))-pwwwww0)));
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(           a*pyyyyy0             );
N2=((n+1)^q2-n^q2)*(    1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))    );
N3=((n+1)^q3-n^q3)*(  1/(e1*Tw)*(-pzzzzz0+e2*pwwwww0-(e*e2*Tw/Ty)*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0)))-pwwwww0 )));
N4=((n+1)^q4-n^q4)*(   1/Ty*(-kp*pwwwww0-ki/a*pxxxxx0-kd*(1/Tab*(pzzzzz0-D*pyyyyy0-(Eq*Vs/xd)*sin(pxxxxx0)-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx0))-pwwwww0)));
for j=1:td/h      
       %                                       f(  x(j) , y(j) , z(j)  delay term replace)
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(   aw*sin(fp*j)+rw*rand(1)+      a*pyyyyy(j)               );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*( aw*sin(fp*j)+rw*rand(1)+     1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))                    );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(  aw*sin(fp*j)+rw*rand(1)+1/(e1*Tw)*(-pzzzzz(j)+e2*pwwwww(j)-(e*e2*Tw/Ty)*(-kp*pwwwww(j)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j))))-pwwwww(j) )) );
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwwwww(j)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))-pwwwww(j))));
       %    f(  x(j) , y(j) , z(j)  delay term replace)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(   aw*sin(fp*j)+rw*rand(1)+     a*pyyyyy(j)                  );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))      );
        N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzzz(j)+e2*pwwwww(j)-(e*e2*Tw/Ty)*(-kp*pwwwww(j)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j))))-pwwwww(j) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwwwww(j)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))-pwwwww(j))));
     
end   
for j=td/h+1:n               
       %                    f(  x(j) , y(j) , z(j) The delay variable is replaced by the value of the previous variable )
       M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    aw*sin(fp*j)+rw*rand(1)+     a*pyyyyy(j)         );
       M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))      );
        M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(aw*sin(fp*j)+rw*rand(1)+  1/(e1*Tw)*(-pzzzzz(j)+e2*pwwwww(j-td/h)-(e*e2*Tw/Ty)*(-kp*pwwwww(j-td/h)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j))))-pwwwww(j-td/h) )));
       M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+ 1/Ty*(-kp*pwwwww(j-td/h)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))-pwwwww(j-td/h))) );
       %                   f(  x(j) , y(j) , z(j)  The delay variable is replaced by the value of the previous variable)
       N1=N1+((n-j+1)^q1-(n-j)^q1)*(   aw*sin(fp*j)+rw*rand(1)+     a*pyyyyy(j)             );
       N2=N2+((n-j+1)^q2-(n-j)^q2)*(  aw*sin(fp*j)+rw*rand(1)+  1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))      );
       N3=N3+((n-j+1)^q3-(n-j)^q3)*(aw*sin(fp*j)+rw*rand(1)+ 1/(e1*Tw)*(-pzzzzz(j)+e2*pwwwww(j-td/h)-(e*e2*Tw/Ty)*(-kp*pwwwww(j-td/h)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j))))-pwwwww(j-td/h) )));
       N4=N4+((n-j+1)^q4-(n-j)^q4)*(aw*sin(fp*j)+rw*rand(1)+  1/Ty*(-kp*pwwwww(j-td/h)-ki/a*pxxxxx(j)-kd*(1/Tab*(pzzzzz(j)-D*pyyyyy(j)-(Eq*Vs/xd)*sin(pxxxxx(j))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(j)))-pwwwww(j-td/h))));
end

pxxxxx1(n+1)=pxxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
pyyyyy1(n+1)=pyyyyy0+h^q2*N2/(gamma(q2)*q2);
pzzzzz1(n+1)=pzzzzz0+h^q3*N3/(gamma(q3)*q3);
pwwwww1(n+1)=pwwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) The delay variable is replaced by the value of the previous variable )
pxxxxx(n+2)=pxxxxx0+h^q1*(        a*pyyyyy(n+1)            +M1)/gamma(q1+2);
pyyyyy(n+2)=pyyyyy0+h^q2*(     1/Tab*(pzzzzz(n+1)-D*pyyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(n+1)))   +M2)/gamma(q2+2);
pzzzzz(n+2)=pzzzzz0+h^q3*(  1/(e1*Tw)*(-pzzzzz(n+1)+e2*pwwwww(n+1-td/h)-(e*e2*Tw/Ty)*(-kp*pwwwww(n+1-td/h)-ki/a*pxxxxx(n+1)-kd*(1/Tab*(pzzzzz(n+1)-D*pyyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(n+1))))-pwwwww(n+1-td/h) ))+M3)/gamma(q3+2);
pwwwww(n+2)=pwwwww0+h^q4*(  1/Ty*(-kp*pwwwww(n+1-td/h)-ki/a*pxxxxx(n+1)-kd*(1/Tab*(pzzzzz(n+1)-D*pyyyyy(n+1)-(Eq*Vs/xd)*sin(pxxxxx(n+1))-(Vs^2/2)*(xd-xq)/(xd*xq)*sin(2*pxxxxx(n+1)))-pwwwww(n+1-td/h)))+M4)/gamma(q4+2);
end

%%
q1=0.88;q2=0.88;q3=0.88;q4=0.88;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mx0=0.01;my0=0.01;mz0=0.01;mw0=0.01;   %initial value
td=0.1;

mx(N+1)=[0];my(N+1)=[0];mz(N+1)=[0];mw(N+1)=[0]; %efficiency need improve
mx1(N+1)=[0];my1(N+1)=[0];mz1(N+1)=[0];mw1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mx1(1)=mx0+h^q1*(   314*my0   -k1*mx0  )/(gamma(q1)*q1);
my1(1)=my0+h^q2*(   1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)   -k2*my0   )/(gamma(q2)*q2);
mz1(1)=mz0+h^q3*(   -2.5*mz0+16.5*mw0   -k3*mz0         )/(gamma(q3)*q3);
mw1(1)=mw0+h^q4*(   -10*mw0   -k4*mw0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mx(1)=mx0+h^q1*((  314*my0-k1*mx0  )+q1*( 314*my0-k1*mx0  ))/gamma(q1+2);
my(1)=my0+h^q2*((  1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0 )+q2*(  1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0  ))/gamma(q2+2);
mz(1)=mz0+h^q3*((  -2.5*mz0+16.5*mw0-k3*mz0  )+q3*(   -2.5*mz0+16.5*mw0-k3*mz0   ))/gamma(q3+2);
mw(1)=mw0+h^q4*(( -10*mw0-k4*mw0  )+q4*(  -10*mw0-k4*mw0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*my0-k1*mx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mz0+16.5*mw0-k3*mz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mw0-k4*mw0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*my0-k1*mx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mz0+16.5*mw0-k3*mz0              );
N4=((n+1)^q4-n^q4)*(      -10*mw0-k4*mw0                     );
mx1(n+1)=mx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
my1(n+1)=my0+h^q2*N2/(gamma(q2)*q2);
mz1(n+1)=mz0+h^q3*N3/(gamma(q3)*q3);
mw1(n+1)=mw0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mx(n+2)=mx0+h^q1*(    314*my(n+1)-k1*mx(n+1)      +M1)/gamma(q1+2);
my(n+2)=my0+h^q2*(     1/19*mz(n+1)-2/19*my(n+1)-1.08/19*sin(mx(n+1))+0.061/19*sin(2*mx(n+1))-k2*my(n+1)      +M2)/gamma(q2+2);
mz(n+2)=mz0+h^q3*(    -2.5*mz(n+1)+16.5*mw(n+1)-k3*mz(n+1)         +M3)/gamma(q3+2);
mw(n+2)=mw0+h^q4*(    -10*mw(n+1)-k4*mw(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*my0-k1*mx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mz0+16.5*mw0-k3*mz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mw0-k4*mw0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*my0-k1*mx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mz0+16.5*mw0-k3*mz0   );
N4=((n+1)^q4-n^q4)*(     -10*mw0-k4*mw0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*my(j)-k1*mx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mz(j)+16.5*mw(j)-k3*mz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mw(j)-k4*mw(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*my(j)-k1*mx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mz(j)+16.5*mw(j)-k3*mz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mw(j)-k4*mw(j)          );
end   
mx1(n+1)=mx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
my1(n+1)=my0+h^q2*N2/(gamma(q2)*q2);
mz1(n+1)=mz0+h^q3*N3/(gamma(q3)*q3);
mw1(n+1)=mw0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mx(n+2)=mx0+h^q1*(       314*my(n+1)-k1*mx(n+1)     +M1)/gamma(q1+2);
my(n+2)=my0+h^q2*(       1/19*mz(n+1)-2/19*my(n+1)-1.08/19*sin(mx(n+1))+0.061/19*sin(2*mx(n+1))-k2*my(n+1)        +M2)/gamma(q2+2);
mz(n+2)=mz0+h^q3*(       -2.5*mz(n+1)+16.5*mw(n+1)-k3*mz(n+1)      +M3)/gamma(q3+2);
mw(n+2)=mw0+h^q4*(     -10*mw(n+1)-k4*mw(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*my0-k1*mx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mz0+16.5*mw0-k3*mz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mw0-k4*mw0       );
N1=((n+1)^q1-n^q1)*(        314*my0-k1*mx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mz0+16.5*mw0-k3*mz0      );
N4=((n+1)^q4-n^q4)*(      -10*mw0-k4*mw0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*my(j)-k1*mx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mz(j)+16.5*mw(j)-k3*mz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mw(j)-k4*mw(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*my(j)-k1*mx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mz(j)+16.5*mw(j)-k3*mz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mw(j)-k4*mw(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*my(j)-k1*mx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mz(j)+16.5*mw(j-td/h)-k3*mz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mw(j-td/h)-k4*mw(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*my(j)-k1*mx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mz(j)+16.5*mw(j-td/h)-k3*mz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mw(j-td/h)-k4*mw(j)          );
end
mx1(n+1)=mx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
my1(n+1)=my0+h^q2*N2/(gamma(q2)*q2);
mz1(n+1)=mz0+h^q3*N3/(gamma(q3)*q3);
mw1(n+1)=mw0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mx(n+2)=mx0+h^q1*(       314*my(n+1)-k1*mx(n+1)      +M1)/gamma(q1+2);
my(n+2)=my0+h^q2*(       1/19*mz(n+1)-2/19*my(n+1)-1.08/19*sin(mx(n+1))+0.061/19*sin(2*mx(n+1))-k2*my(n+1)       +M2)/gamma(q2+2);
mz(n+2)=mz0+h^q3*(       -2.5*mz(n+1)+16.5*mw(n+1-td/h)-k3*mz(n+1)    +M3)/gamma(q3+2);
mw(n+2)=mw0+h^q4*(      -10*mw(n+1-td/h)-k4*mw(n+1)    +M4)/gamma(q4+2);
end


q1=0.9;q2=0.9;q3=0.9;q4=0.9;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mxx0=0.01;myy0=0.01;mzz0=0.01;mww0=0.01;   %initial value
td=0.1;

mxx(N+1)=[0];myy(N+1)=[0];mzz(N+1)=[0];mww(N+1)=[0]; %efficiency need improve
mxx1(N+1)=[0];myy1(N+1)=[0];mzz1(N+1)=[0];mww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mxx1(1)=mxx0+h^q1*(   314*myy0   -k1*mxx0  )/(gamma(q1)*q1);
myy1(1)=myy0+h^q2*(   1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)   -k2*myy0   )/(gamma(q2)*q2);
mzz1(1)=mzz0+h^q3*(   -2.5*mzz0+16.5*mww0   -k3*mzz0         )/(gamma(q3)*q3);
mww1(1)=mww0+h^q4*(   -10*mww0   -k4*mww0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mxx(1)=mxx0+h^q1*((  314*myy0-k1*mxx0  )+q1*( 314*myy0-k1*mxx0  ))/gamma(q1+2);
myy(1)=myy0+h^q2*((  1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0 )+q2*(  1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0  ))/gamma(q2+2);
mzz(1)=mzz0+h^q3*((  -2.5*mzz0+16.5*mww0-k3*mzz0  )+q3*(   -2.5*mzz0+16.5*mww0-k3*mzz0   ))/gamma(q3+2);
mww(1)=mww0+h^q4*(( -10*mww0-k4*mww0  )+q4*(  -10*mww0-k4*mww0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myy0-k1*mxx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mzz0+16.5*mww0-k3*mzz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mww0-k4*mww0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*myy0-k1*mxx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mzz0+16.5*mww0-k3*mzz0              );
N4=((n+1)^q4-n^q4)*(      -10*mww0-k4*mww0                     );
mxx1(n+1)=mxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myy1(n+1)=myy0+h^q2*N2/(gamma(q2)*q2);
mzz1(n+1)=mzz0+h^q3*N3/(gamma(q3)*q3);
mww1(n+1)=mww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mxx(n+2)=mxx0+h^q1*(    314*myy(n+1)-k1*mxx(n+1)      +M1)/gamma(q1+2);
myy(n+2)=myy0+h^q2*(     1/19*mzz(n+1)-2/19*myy(n+1)-1.08/19*sin(mxx(n+1))+0.061/19*sin(2*mxx(n+1))-k2*myy(n+1)      +M2)/gamma(q2+2);
mzz(n+2)=mzz0+h^q3*(    -2.5*mzz(n+1)+16.5*mww(n+1)-k3*mzz(n+1)         +M3)/gamma(q3+2);
mww(n+2)=mww0+h^q4*(    -10*mww(n+1)-k4*mww(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myy0-k1*mxx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mzz0+16.5*mww0-k3*mzz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mww0-k4*mww0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*myy0-k1*mxx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mzz0+16.5*mww0-k3*mzz0   );
N4=((n+1)^q4-n^q4)*(     -10*mww0-k4*mww0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*myy(j)-k1*mxx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mzz(j)+16.5*mww(j)-k3*mzz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mww(j)-k4*mww(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myy(j)-k1*mxx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzz(j)+16.5*mww(j)-k3*mzz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mww(j)-k4*mww(j)          );
end   
mxx1(n+1)=mxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myy1(n+1)=myy0+h^q2*N2/(gamma(q2)*q2);
mzz1(n+1)=mzz0+h^q3*N3/(gamma(q3)*q3);
mww1(n+1)=mww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxx(n+2)=mxx0+h^q1*(       314*myy(n+1)-k1*mxx(n+1)     +M1)/gamma(q1+2);
myy(n+2)=myy0+h^q2*(       1/19*mzz(n+1)-2/19*myy(n+1)-1.08/19*sin(mxx(n+1))+0.061/19*sin(2*mxx(n+1))-k2*myy(n+1)        +M2)/gamma(q2+2);
mzz(n+2)=mzz0+h^q3*(       -2.5*mzz(n+1)+16.5*mww(n+1)-k3*mzz(n+1)      +M3)/gamma(q3+2);
mww(n+2)=mww0+h^q4*(     -10*mww(n+1)-k4*mww(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*myy0-k1*mxx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mzz0+16.5*mww0-k3*mzz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mww0-k4*mww0       );
N1=((n+1)^q1-n^q1)*(        314*myy0-k1*mxx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mzz0+16.5*mww0-k3*mzz0      );
N4=((n+1)^q4-n^q4)*(      -10*mww0-k4*mww0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*myy(j)-k1*mxx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mzz(j)+16.5*mww(j)-k3*mzz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mww(j)-k4*mww(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*myy(j)-k1*mxx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzz(j)+16.5*mww(j)-k3*mzz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mww(j)-k4*mww(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*myy(j)-k1*mxx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mzz(j)+16.5*mww(j-td/h)-k3*mzz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mww(j-td/h)-k4*mww(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myy(j)-k1*mxx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mzz(j)+16.5*mww(j-td/h)-k3*mzz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mww(j-td/h)-k4*mww(j)          );
end
mxx1(n+1)=mxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myy1(n+1)=myy0+h^q2*N2/(gamma(q2)*q2);
mzz1(n+1)=mzz0+h^q3*N3/(gamma(q3)*q3);
mww1(n+1)=mww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxx(n+2)=mxx0+h^q1*(       314*myy(n+1)-k1*mxx(n+1)      +M1)/gamma(q1+2);
myy(n+2)=myy0+h^q2*(       1/19*mzz(n+1)-2/19*myy(n+1)-1.08/19*sin(mxx(n+1))+0.061/19*sin(2*mxx(n+1))-k2*myy(n+1)       +M2)/gamma(q2+2);
mzz(n+2)=mzz0+h^q3*(       -2.5*mzz(n+1)+16.5*mww(n+1-td/h)-k3*mzz(n+1)    +M3)/gamma(q3+2);
mww(n+2)=mww0+h^q4*(      -10*mww(n+1-td/h)-k4*mww(n+1)    +M4)/gamma(q4+2);
end


q1=0.92;q2=0.92;q3=0.92;q4=0.92;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mxxx0=0.01;myyy0=0.01;mzzz0=0.01;mwww0=0.01;   %initial value
td=0.1;

mxxx(N+1)=[0];myyy(N+1)=[0];mzzz(N+1)=[0];mwww(N+1)=[0]; %efficiency need improve
mxxx1(N+1)=[0];myyy1(N+1)=[0];mzzz1(N+1)=[0];mwww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mxxx1(1)=mxxx0+h^q1*(   314*myyy0   -k1*mxxx0  )/(gamma(q1)*q1);
myyy1(1)=myyy0+h^q2*(   1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)   -k2*myyy0   )/(gamma(q2)*q2);
mzzz1(1)=mzzz0+h^q3*(   -2.5*mzzz0+16.5*mwww0   -k3*mzzz0         )/(gamma(q3)*q3);
mwww1(1)=mwww0+h^q4*(   -10*mwww0   -k4*mwww0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mxxx(1)=mxxx0+h^q1*((  314*myyy0-k1*mxxx0  )+q1*( 314*myyy0-k1*mxxx0  ))/gamma(q1+2);
myyy(1)=myyy0+h^q2*((  1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0 )+q2*(  1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0  ))/gamma(q2+2);
mzzz(1)=mzzz0+h^q3*((  -2.5*mzzz0+16.5*mwww0-k3*mzzz0  )+q3*(   -2.5*mzzz0+16.5*mwww0-k3*mzzz0   ))/gamma(q3+2);
mwww(1)=mwww0+h^q4*(( -10*mwww0-k4*mwww0  )+q4*(  -10*mwww0-k4*mwww0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyy0-k1*mxxx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mzzz0+16.5*mwww0-k3*mzzz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mwww0-k4*mwww0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*myyy0-k1*mxxx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mzzz0+16.5*mwww0-k3*mzzz0              );
N4=((n+1)^q4-n^q4)*(      -10*mwww0-k4*mwww0                     );
mxxx1(n+1)=mxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyy1(n+1)=myyy0+h^q2*N2/(gamma(q2)*q2);
mzzz1(n+1)=mzzz0+h^q3*N3/(gamma(q3)*q3);
mwww1(n+1)=mwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mxxx(n+2)=mxxx0+h^q1*(    314*myyy(n+1)-k1*mxxx(n+1)      +M1)/gamma(q1+2);
myyy(n+2)=myyy0+h^q2*(     1/19*mzzz(n+1)-2/19*myyy(n+1)-1.08/19*sin(mxxx(n+1))+0.061/19*sin(2*mxxx(n+1))-k2*myyy(n+1)      +M2)/gamma(q2+2);
mzzz(n+2)=mzzz0+h^q3*(    -2.5*mzzz(n+1)+16.5*mwww(n+1)-k3*mzzz(n+1)         +M3)/gamma(q3+2);
mwww(n+2)=mwww0+h^q4*(    -10*mwww(n+1)-k4*mwww(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyy0-k1*mxxx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mzzz0+16.5*mwww0-k3*mzzz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mwww0-k4*mwww0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*myyy0-k1*mxxx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mzzz0+16.5*mwww0-k3*mzzz0   );
N4=((n+1)^q4-n^q4)*(     -10*mwww0-k4*mwww0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*myyy(j)-k1*mxxx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mzzz(j)+16.5*mwww(j)-k3*mzzz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mwww(j)-k4*mwww(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyy(j)-k1*mxxx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzz(j)+16.5*mwww(j)-k3*mzzz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mwww(j)-k4*mwww(j)          );
end   
mxxx1(n+1)=mxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyy1(n+1)=myyy0+h^q2*N2/(gamma(q2)*q2);
mzzz1(n+1)=mzzz0+h^q3*N3/(gamma(q3)*q3);
mwww1(n+1)=mwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxx(n+2)=mxxx0+h^q1*(       314*myyy(n+1)-k1*mxxx(n+1)     +M1)/gamma(q1+2);
myyy(n+2)=myyy0+h^q2*(       1/19*mzzz(n+1)-2/19*myyy(n+1)-1.08/19*sin(mxxx(n+1))+0.061/19*sin(2*mxxx(n+1))-k2*myyy(n+1)        +M2)/gamma(q2+2);
mzzz(n+2)=mzzz0+h^q3*(       -2.5*mzzz(n+1)+16.5*mwww(n+1)-k3*mzzz(n+1)      +M3)/gamma(q3+2);
mwww(n+2)=mwww0+h^q4*(     -10*mwww(n+1)-k4*mwww(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*myyy0-k1*mxxx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mzzz0+16.5*mwww0-k3*mzzz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mwww0-k4*mwww0       );
N1=((n+1)^q1-n^q1)*(        314*myyy0-k1*mxxx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mzzz0+16.5*mwww0-k3*mzzz0      );
N4=((n+1)^q4-n^q4)*(      -10*mwww0-k4*mwww0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*myyy(j)-k1*mxxx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mzzz(j)+16.5*mwww(j)-k3*mzzz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mwww(j)-k4*mwww(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*myyy(j)-k1*mxxx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzz(j)+16.5*mwww(j)-k3*mzzz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwww(j)-k4*mwww(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*myyy(j)-k1*mxxx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mzzz(j)+16.5*mwww(j-td/h)-k3*mzzz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mwww(j-td/h)-k4*mwww(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyy(j)-k1*mxxx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mzzz(j)+16.5*mwww(j-td/h)-k3*mzzz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwww(j-td/h)-k4*mwww(j)          );
end
mxxx1(n+1)=mxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyy1(n+1)=myyy0+h^q2*N2/(gamma(q2)*q2);
mzzz1(n+1)=mzzz0+h^q3*N3/(gamma(q3)*q3);
mwww1(n+1)=mwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxx(n+2)=mxxx0+h^q1*(       314*myyy(n+1)-k1*mxxx(n+1)      +M1)/gamma(q1+2);
myyy(n+2)=myyy0+h^q2*(       1/19*mzzz(n+1)-2/19*myyy(n+1)-1.08/19*sin(mxxx(n+1))+0.061/19*sin(2*mxxx(n+1))-k2*myyy(n+1)       +M2)/gamma(q2+2);
mzzz(n+2)=mzzz0+h^q3*(       -2.5*mzzz(n+1)+16.5*mwww(n+1-td/h)-k3*mzzz(n+1)    +M3)/gamma(q3+2);
mwww(n+2)=mwww0+h^q4*(      -10*mwww(n+1-td/h)-k4*mwww(n+1)    +M4)/gamma(q4+2);
end


q1=0.94;q2=0.94;q3=0.94;q4=0.94;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mxxxx0=0.01;myyyy0=0.01;mzzzz0=0.01;mwwww0=0.01;   %initial value
td=0.1;

mxxxx(N+1)=[0];myyyy(N+1)=[0];mzzzz(N+1)=[0];mwwww(N+1)=[0]; %efficiency need improve
mxxxx1(N+1)=[0];myyyy1(N+1)=[0];mzzzz1(N+1)=[0];mwwww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mxxxx1(1)=mxxxx0+h^q1*(   314*myyyy0   -k1*mxxxx0  )/(gamma(q1)*q1);
myyyy1(1)=myyyy0+h^q2*(   1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)   -k2*myyyy0   )/(gamma(q2)*q2);
mzzzz1(1)=mzzzz0+h^q3*(   -2.5*mzzzz0+16.5*mwwww0   -k3*mzzzz0         )/(gamma(q3)*q3);
mwwww1(1)=mwwww0+h^q4*(   -10*mwwww0   -k4*mwwww0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mxxxx(1)=mxxxx0+h^q1*((  314*myyyy0-k1*mxxxx0  )+q1*( 314*myyyy0-k1*mxxxx0  ))/gamma(q1+2);
myyyy(1)=myyyy0+h^q2*((  1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0 )+q2*(  1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0  ))/gamma(q2+2);
mzzzz(1)=mzzzz0+h^q3*((  -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0  )+q3*(   -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0   ))/gamma(q3+2);
mwwww(1)=mwwww0+h^q4*(( -10*mwwww0-k4*mwwww0  )+q4*(  -10*mwwww0-k4*mwwww0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyyy0-k1*mxxxx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mwwww0-k4*mwwww0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*myyyy0-k1*mxxxx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0              );
N4=((n+1)^q4-n^q4)*(      -10*mwwww0-k4*mwwww0                     );
mxxxx1(n+1)=mxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyy1(n+1)=myyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzz1(n+1)=mzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwww1(n+1)=mwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mxxxx(n+2)=mxxxx0+h^q1*(    314*myyyy(n+1)-k1*mxxxx(n+1)      +M1)/gamma(q1+2);
myyyy(n+2)=myyyy0+h^q2*(     1/19*mzzzz(n+1)-2/19*myyyy(n+1)-1.08/19*sin(mxxxx(n+1))+0.061/19*sin(2*mxxxx(n+1))-k2*myyyy(n+1)      +M2)/gamma(q2+2);
mzzzz(n+2)=mzzzz0+h^q3*(    -2.5*mzzzz(n+1)+16.5*mwwww(n+1)-k3*mzzzz(n+1)         +M3)/gamma(q3+2);
mwwww(n+2)=mwwww0+h^q4*(    -10*mwwww(n+1)-k4*mwwww(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyyy0-k1*mxxxx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mwwww0-k4*mwwww0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*myyyy0-k1*mxxxx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0   );
N4=((n+1)^q4-n^q4)*(     -10*mwwww0-k4*mwwww0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*myyyy(j)-k1*mxxxx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mzzzz(j)+16.5*mwwww(j)-k3*mzzzz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mwwww(j)-k4*mwwww(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyyy(j)-k1*mxxxx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzzz(j)+16.5*mwwww(j)-k3*mzzzz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mwwww(j)-k4*mwwww(j)          );
end   
mxxxx1(n+1)=mxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyy1(n+1)=myyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzz1(n+1)=mzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwww1(n+1)=mwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxxx(n+2)=mxxxx0+h^q1*(       314*myyyy(n+1)-k1*mxxxx(n+1)     +M1)/gamma(q1+2);
myyyy(n+2)=myyyy0+h^q2*(       1/19*mzzzz(n+1)-2/19*myyyy(n+1)-1.08/19*sin(mxxxx(n+1))+0.061/19*sin(2*mxxxx(n+1))-k2*myyyy(n+1)        +M2)/gamma(q2+2);
mzzzz(n+2)=mzzzz0+h^q3*(       -2.5*mzzzz(n+1)+16.5*mwwww(n+1)-k3*mzzzz(n+1)      +M3)/gamma(q3+2);
mwwww(n+2)=mwwww0+h^q4*(     -10*mwwww(n+1)-k4*mwwww(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*myyyy0-k1*mxxxx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mwwww0-k4*mwwww0       );
N1=((n+1)^q1-n^q1)*(        314*myyyy0-k1*mxxxx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0      );
N4=((n+1)^q4-n^q4)*(      -10*mwwww0-k4*mwwww0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*myyyy(j)-k1*mxxxx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mzzzz(j)+16.5*mwwww(j)-k3*mzzzz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mwwww(j)-k4*mwwww(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*myyyy(j)-k1*mxxxx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzzz(j)+16.5*mwwww(j)-k3*mzzzz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwwww(j)-k4*mwwww(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*myyyy(j)-k1*mxxxx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mzzzz(j)+16.5*mwwww(j-td/h)-k3*mzzzz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mwwww(j-td/h)-k4*mwwww(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyyy(j)-k1*mxxxx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mzzzz(j)+16.5*mwwww(j-td/h)-k3*mzzzz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwwww(j-td/h)-k4*mwwww(j)          );
end
mxxxx1(n+1)=mxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyy1(n+1)=myyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzz1(n+1)=mzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwww1(n+1)=mwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxxx(n+2)=mxxxx0+h^q1*(       314*myyyy(n+1)-k1*mxxxx(n+1)      +M1)/gamma(q1+2);
myyyy(n+2)=myyyy0+h^q2*(       1/19*mzzzz(n+1)-2/19*myyyy(n+1)-1.08/19*sin(mxxxx(n+1))+0.061/19*sin(2*mxxxx(n+1))-k2*myyyy(n+1)       +M2)/gamma(q2+2);
mzzzz(n+2)=mzzzz0+h^q3*(       -2.5*mzzzz(n+1)+16.5*mwwww(n+1-td/h)-k3*mzzzz(n+1)    +M3)/gamma(q3+2);
mwwww(n+2)=mwwww0+h^q4*(      -10*mwwww(n+1-td/h)-k4*mwwww(n+1)    +M4)/gamma(q4+2);
end


q1=0.95;q2=0.95;q3=0.95;q4=0.95;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mxxxxx0=0.01;myyyyy0=0.01;mzzzzz0=0.01;mwwwww0=0.01;   %initial value
td=0.08;

mxxxxx(N+1)=[0];myyyyy(N+1)=[0];mzzzzz(N+1)=[0];mwwwww(N+1)=[0]; %efficiency need improve
mxxxxx1(N+1)=[0];myyyyy1(N+1)=[0];mzzzzz1(N+1)=[0];mwwwww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mxxxxx1(1)=mxxxxx0+h^q1*(   314*myyyyy0   -k1*mxxxxx0  )/(gamma(q1)*q1);
myyyyy1(1)=myyyyy0+h^q2*(   1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)   -k2*myyyyy0   )/(gamma(q2)*q2);
mzzzzz1(1)=mzzzzz0+h^q3*(   -2.5*mzzzzz0+16.5*mwwwww0   -k3*mzzzzz0         )/(gamma(q3)*q3);
mwwwww1(1)=mwwwww0+h^q4*(   -10*mwwwww0   -k4*mwwwww0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mxxxxx(1)=mxxxxx0+h^q1*((  314*myyyyy0-k1*mxxxxx0  )+q1*( 314*myyyyy0-k1*mxxxxx0  ))/gamma(q1+2);
myyyyy(1)=myyyyy0+h^q2*((  1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0 )+q2*(  1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0  ))/gamma(q2+2);
mzzzzz(1)=mzzzzz0+h^q3*((  -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0  )+q3*(   -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0   ))/gamma(q3+2);
mwwwww(1)=mwwwww0+h^q4*(( -10*mwwwww0-k4*mwwwww0  )+q4*(  -10*mwwwww0-k4*mwwwww0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyyyy0-k1*mxxxxx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mwwwww0-k4*mwwwww0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*myyyyy0-k1*mxxxxx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0              );
N4=((n+1)^q4-n^q4)*(      -10*mwwwww0-k4*mwwwww0                     );
mxxxxx1(n+1)=mxxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyyy1(n+1)=myyyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzzz1(n+1)=mzzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwwww1(n+1)=mwwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mxxxxx(n+2)=mxxxxx0+h^q1*(    314*myyyyy(n+1)-k1*mxxxxx(n+1)      +M1)/gamma(q1+2);
myyyyy(n+2)=myyyyy0+h^q2*(     1/19*mzzzzz(n+1)-2/19*myyyyy(n+1)-1.08/19*sin(mxxxxx(n+1))+0.061/19*sin(2*mxxxxx(n+1))-k2*myyyyy(n+1)      +M2)/gamma(q2+2);
mzzzzz(n+2)=mzzzzz0+h^q3*(    -2.5*mzzzzz(n+1)+16.5*mwwwww(n+1)-k3*mzzzzz(n+1)         +M3)/gamma(q3+2);
mwwwww(n+2)=mwwwww0+h^q4*(    -10*mwwwww(n+1)-k4*mwwwww(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyyyy0-k1*mxxxxx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mwwwww0-k4*mwwwww0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*myyyyy0-k1*mxxxxx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0   );
N4=((n+1)^q4-n^q4)*(     -10*mwwwww0-k4*mwwwww0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*myyyyy(j)-k1*mxxxxx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mzzzzz(j)+16.5*mwwwww(j)-k3*mzzzzz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mwwwww(j)-k4*mwwwww(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyyyy(j)-k1*mxxxxx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzzzz(j)+16.5*mwwwww(j)-k3*mzzzzz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mwwwww(j)-k4*mwwwww(j)          );
end   
mxxxxx1(n+1)=mxxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyyy1(n+1)=myyyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzzz1(n+1)=mzzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwwww1(n+1)=mwwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxxxx(n+2)=mxxxxx0+h^q1*(       314*myyyyy(n+1)-k1*mxxxxx(n+1)     +M1)/gamma(q1+2);
myyyyy(n+2)=myyyyy0+h^q2*(       1/19*mzzzzz(n+1)-2/19*myyyyy(n+1)-1.08/19*sin(mxxxxx(n+1))+0.061/19*sin(2*mxxxxx(n+1))-k2*myyyyy(n+1)        +M2)/gamma(q2+2);
mzzzzz(n+2)=mzzzzz0+h^q3*(       -2.5*mzzzzz(n+1)+16.5*mwwwww(n+1)-k3*mzzzzz(n+1)      +M3)/gamma(q3+2);
mwwwww(n+2)=mwwwww0+h^q4*(     -10*mwwwww(n+1)-k4*mwwwww(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*myyyyy0-k1*mxxxxx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mwwwww0-k4*mwwwww0       );
N1=((n+1)^q1-n^q1)*(        314*myyyyy0-k1*mxxxxx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0      );
N4=((n+1)^q4-n^q4)*(      -10*mwwwww0-k4*mwwwww0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*myyyyy(j)-k1*mxxxxx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mzzzzz(j)+16.5*mwwwww(j)-k3*mzzzzz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mwwwww(j)-k4*mwwwww(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*myyyyy(j)-k1*mxxxxx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzzzz(j)+16.5*mwwwww(j)-k3*mzzzzz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwwwww(j)-k4*mwwwww(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*myyyyy(j)-k1*mxxxxx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mzzzzz(j)+16.5*mwwwww(j-td/h)-k3*mzzzzz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mwwwww(j-td/h)-k4*mwwwww(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyyyy(j)-k1*mxxxxx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mzzzzz(j)+16.5*mwwwww(j-td/h)-k3*mzzzzz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwwwww(j-td/h)-k4*mwwwww(j)          );
end
mxxxxx1(n+1)=mxxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyyy1(n+1)=myyyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzzz1(n+1)=mzzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwwww1(n+1)=mwwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxxxx(n+2)=mxxxxx0+h^q1*(       314*myyyyy(n+1)-k1*mxxxxx(n+1)      +M1)/gamma(q1+2);
myyyyy(n+2)=myyyyy0+h^q2*(       1/19*mzzzzz(n+1)-2/19*myyyyy(n+1)-1.08/19*sin(mxxxxx(n+1))+0.061/19*sin(2*mxxxxx(n+1))-k2*myyyyy(n+1)       +M2)/gamma(q2+2);
mzzzzz(n+2)=mzzzzz0+h^q3*(       -2.5*mzzzzz(n+1)+16.5*mwwwww(n+1-td/h)-k3*mzzzzz(n+1)    +M3)/gamma(q3+2);
mwwwww(n+2)=mwwwww0+h^q4*(      -10*mwwwww(n+1-td/h)-k4*mwwwww(n+1)    +M4)/gamma(q4+2);
end
%%
A11=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 17231/16951 -2/9 0;0 0 0 0 0 0];
A12=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 1577/16951 -2/9 0;0 0 0 0 0 0];
Ad=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 -7/45;0 0 0 0 0 -10];B=[0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];I=eye(6);
setlmis([]);yip1=0.1;yip2=0.1;dert1=0.1;dert2=0.1;Q=lmivar(1,[6 1]);M1=lmivar(2,[3 6]);M2=lmivar(2,[3 6]);lmiterm([1 1 1 Q],A11,1,'s');lmiterm([1 1 1 M1],B,1,'s');lmiterm([1 1 1 Q],1,1);
lmiterm([1 1 1 0],yip1*(Ad*Ad'));lmiterm([1 1 2 0],0);lmiterm([1 1 3 0],0);lmiterm([1 2 2 Q],-1,1);lmiterm([1 2 3 Q],1,1);lmiterm([1 3 3 0],-yip1*I);lmiterm([2 1 1 Q],A12,1,'s');
lmiterm([2 1 1 M2],B,1,'s');lmiterm([2 1 1 Q],1,1);lmiterm([2 1 1 0],yip2*(Ad*Ad'));lmiterm([2 1 2 0],0);lmiterm([2 1 3 0],0);lmiterm([2 2 2 Q],-1,1);lmiterm([2 2 3 Q],1,1);lmiterm([2 3 3 0],-yip2*I);
lmiterm([3 1 1 Q],1,A11'+A12','s');lmiterm([3 1 1 -M1],1,B','s');lmiterm([3 1 1 -M2],1,B','s');lmiterm([3 1 1 Q],2,1);lmiterm([3 1 1 0],(dert1+dert2)*(Ad*Ad'));lmiterm([3 1 2 0],0);
lmiterm([3 1 3 0],0);lmiterm([3 1 4 0],0);lmiterm([3 2 2 Q],-2,1);lmiterm([3 2 3 Q],1,1);lmiterm([3 2 4 Q],1,1);lmiterm([3 3 3 0],-dert1*I);lmiterm([3 3 4 0],0);lmiterm([3 4 4 0],-dert2*I);
lmiterm([-4 1 1 Q],1,1);lmisys=getlmis;[tmin,xfeas]=feasp(lmisys);MATQ=dec2mat(lmisys,xfeas,Q);MATM1=dec2mat(lmisys,xfeas,M1);MATM2=dec2mat(lmisys,xfeas,M2);k1=MATM1/MATQ;k2=MATM2/MATQ;
q1=0.88;q2=0.88;q3=0.88;q4=0.88;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mx0=0.01;my0=0.01;mz0=0.01;mw0=0.01;   %initial value
td=0.1;

mx(N+1)=[0];my(N+1)=[0];mz(N+1)=[0];mw(N+1)=[0]; %efficiency need improve
mx1(N+1)=[0];my1(N+1)=[0];mz1(N+1)=[0];mw1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mx1(1)=mx0+h^q1*(   314*my0   -k1*mx0  )/(gamma(q1)*q1);
my1(1)=my0+h^q2*(   1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)   -k2*my0   )/(gamma(q2)*q2);
mz1(1)=mz0+h^q3*(   -2.5*mz0+16.5*mw0   -k3*mz0         )/(gamma(q3)*q3);
mw1(1)=mw0+h^q4*(   -10*mw0   -k4*mw0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mx(1)=mx0+h^q1*((  314*my0-k1*mx0  )+q1*( 314*my0-k1*mx0  ))/gamma(q1+2);
my(1)=my0+h^q2*((  1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0 )+q2*(  1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0  ))/gamma(q2+2);
mz(1)=mz0+h^q3*((  -2.5*mz0+16.5*mw0-k3*mz0  )+q3*(   -2.5*mz0+16.5*mw0-k3*mz0   ))/gamma(q3+2);
mw(1)=mw0+h^q4*(( -10*mw0-k4*mw0  )+q4*(  -10*mw0-k4*mw0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*my0-k1*mx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mz0+16.5*mw0-k3*mz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mw0-k4*mw0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*my0-k1*mx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mz0+16.5*mw0-k3*mz0              );
N4=((n+1)^q4-n^q4)*(      -10*mw0-k4*mw0                     );
mx1(n+1)=mx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
my1(n+1)=my0+h^q2*N2/(gamma(q2)*q2);
mz1(n+1)=mz0+h^q3*N3/(gamma(q3)*q3);
mw1(n+1)=mw0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mx(n+2)=mx0+h^q1*(    314*my(n+1)-k1*mx(n+1)      +M1)/gamma(q1+2);
my(n+2)=my0+h^q2*(     1/19*mz(n+1)-2/19*my(n+1)-1.08/19*sin(mx(n+1))+0.061/19*sin(2*mx(n+1))-k2*my(n+1)      +M2)/gamma(q2+2);
mz(n+2)=mz0+h^q3*(    -2.5*mz(n+1)+16.5*mw(n+1)-k3*mz(n+1)         +M3)/gamma(q3+2);
mw(n+2)=mw0+h^q4*(    -10*mw(n+1)-k4*mw(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*my0-k1*mx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mz0+16.5*mw0-k3*mz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mw0-k4*mw0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*my0-k1*mx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mz0+16.5*mw0-k3*mz0   );
N4=((n+1)^q4-n^q4)*(     -10*mw0-k4*mw0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*my(j)-k1*mx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mz(j)+16.5*mw(j)-k3*mz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mw(j)-k4*mw(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*my(j)-k1*mx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mz(j)+16.5*mw(j)-k3*mz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mw(j)-k4*mw(j)          );
end   
mx1(n+1)=mx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
my1(n+1)=my0+h^q2*N2/(gamma(q2)*q2);
mz1(n+1)=mz0+h^q3*N3/(gamma(q3)*q3);
mw1(n+1)=mw0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mx(n+2)=mx0+h^q1*(       314*my(n+1)-k1*mx(n+1)     +M1)/gamma(q1+2);
my(n+2)=my0+h^q2*(       1/19*mz(n+1)-2/19*my(n+1)-1.08/19*sin(mx(n+1))+0.061/19*sin(2*mx(n+1))-k2*my(n+1)        +M2)/gamma(q2+2);
mz(n+2)=mz0+h^q3*(       -2.5*mz(n+1)+16.5*mw(n+1)-k3*mz(n+1)      +M3)/gamma(q3+2);
mw(n+2)=mw0+h^q4*(     -10*mw(n+1)-k4*mw(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*my0-k1*mx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mz0+16.5*mw0-k3*mz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mw0-k4*mw0       );
N1=((n+1)^q1-n^q1)*(        314*my0-k1*mx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mz0-2/19*my0-1.08/19*sin(mx0)+0.061/19*sin(2*mx0)-k2*my0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mz0+16.5*mw0-k3*mz0      );
N4=((n+1)^q4-n^q4)*(      -10*mw0-k4*mw0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*my(j)-k1*mx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mz(j)+16.5*mw(j)-k3*mz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mw(j)-k4*mw(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*my(j)-k1*mx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mz(j)+16.5*mw(j)-k3*mz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mw(j)-k4*mw(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*my(j)-k1*mx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mz(j)+16.5*mw(j-td/h)-k3*mz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mw(j-td/h)-k4*mw(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*my(j)-k1*mx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mz(j)-2/19*my(j)-1.08/19*sin(mx(j))+0.061/19*sin(2*mx(j))-k2*my(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mz(j)+16.5*mw(j-td/h)-k3*mz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mw(j-td/h)-k4*mw(j)          );
end
mx1(n+1)=mx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
my1(n+1)=my0+h^q2*N2/(gamma(q2)*q2);
mz1(n+1)=mz0+h^q3*N3/(gamma(q3)*q3);
mw1(n+1)=mw0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mx(n+2)=mx0+h^q1*(       314*my(n+1)-k1*mx(n+1)      +M1)/gamma(q1+2);
my(n+2)=my0+h^q2*(       1/19*mz(n+1)-2/19*my(n+1)-1.08/19*sin(mx(n+1))+0.061/19*sin(2*mx(n+1))-k2*my(n+1)       +M2)/gamma(q2+2);
mz(n+2)=mz0+h^q3*(       -2.5*mz(n+1)+16.5*mw(n+1-td/h)-k3*mz(n+1)    +M3)/gamma(q3+2);
mw(n+2)=mw0+h^q4*(      -10*mw(n+1-td/h)-k4*mw(n+1)    +M4)/gamma(q4+2);
end


q1=0.9;q2=0.9;q3=0.9;q4=0.9;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mxx0=0.01;myy0=0.01;mzz0=0.01;mww0=0.01;   %initial value
td=0.1;

mxx(N+1)=[0];myy(N+1)=[0];mzz(N+1)=[0];mww(N+1)=[0]; %efficiency need improve
mxx1(N+1)=[0];myy1(N+1)=[0];mzz1(N+1)=[0];mww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mxx1(1)=mxx0+h^q1*(   314*myy0   -k1*mxx0  )/(gamma(q1)*q1);
myy1(1)=myy0+h^q2*(   1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)   -k2*myy0   )/(gamma(q2)*q2);
mzz1(1)=mzz0+h^q3*(   -2.5*mzz0+16.5*mww0   -k3*mzz0         )/(gamma(q3)*q3);
mww1(1)=mww0+h^q4*(   -10*mww0   -k4*mww0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mxx(1)=mxx0+h^q1*((  314*myy0-k1*mxx0  )+q1*( 314*myy0-k1*mxx0  ))/gamma(q1+2);
myy(1)=myy0+h^q2*((  1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0 )+q2*(  1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0  ))/gamma(q2+2);
mzz(1)=mzz0+h^q3*((  -2.5*mzz0+16.5*mww0-k3*mzz0  )+q3*(   -2.5*mzz0+16.5*mww0-k3*mzz0   ))/gamma(q3+2);
mww(1)=mww0+h^q4*(( -10*mww0-k4*mww0  )+q4*(  -10*mww0-k4*mww0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myy0-k1*mxx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mzz0+16.5*mww0-k3*mzz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mww0-k4*mww0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*myy0-k1*mxx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mzz0+16.5*mww0-k3*mzz0              );
N4=((n+1)^q4-n^q4)*(      -10*mww0-k4*mww0                     );
mxx1(n+1)=mxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myy1(n+1)=myy0+h^q2*N2/(gamma(q2)*q2);
mzz1(n+1)=mzz0+h^q3*N3/(gamma(q3)*q3);
mww1(n+1)=mww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mxx(n+2)=mxx0+h^q1*(    314*myy(n+1)-k1*mxx(n+1)      +M1)/gamma(q1+2);
myy(n+2)=myy0+h^q2*(     1/19*mzz(n+1)-2/19*myy(n+1)-1.08/19*sin(mxx(n+1))+0.061/19*sin(2*mxx(n+1))-k2*myy(n+1)      +M2)/gamma(q2+2);
mzz(n+2)=mzz0+h^q3*(    -2.5*mzz(n+1)+16.5*mww(n+1)-k3*mzz(n+1)         +M3)/gamma(q3+2);
mww(n+2)=mww0+h^q4*(    -10*mww(n+1)-k4*mww(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myy0-k1*mxx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mzz0+16.5*mww0-k3*mzz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mww0-k4*mww0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*myy0-k1*mxx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mzz0+16.5*mww0-k3*mzz0   );
N4=((n+1)^q4-n^q4)*(     -10*mww0-k4*mww0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*myy(j)-k1*mxx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mzz(j)+16.5*mww(j)-k3*mzz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mww(j)-k4*mww(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myy(j)-k1*mxx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzz(j)+16.5*mww(j)-k3*mzz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mww(j)-k4*mww(j)          );
end   
mxx1(n+1)=mxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myy1(n+1)=myy0+h^q2*N2/(gamma(q2)*q2);
mzz1(n+1)=mzz0+h^q3*N3/(gamma(q3)*q3);
mww1(n+1)=mww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxx(n+2)=mxx0+h^q1*(       314*myy(n+1)-k1*mxx(n+1)     +M1)/gamma(q1+2);
myy(n+2)=myy0+h^q2*(       1/19*mzz(n+1)-2/19*myy(n+1)-1.08/19*sin(mxx(n+1))+0.061/19*sin(2*mxx(n+1))-k2*myy(n+1)        +M2)/gamma(q2+2);
mzz(n+2)=mzz0+h^q3*(       -2.5*mzz(n+1)+16.5*mww(n+1)-k3*mzz(n+1)      +M3)/gamma(q3+2);
mww(n+2)=mww0+h^q4*(     -10*mww(n+1)-k4*mww(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*myy0-k1*mxx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mzz0+16.5*mww0-k3*mzz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mww0-k4*mww0       );
N1=((n+1)^q1-n^q1)*(        314*myy0-k1*mxx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mzz0-2/19*myy0-1.08/19*sin(mxx0)+0.061/19*sin(2*mxx0)-k2*myy0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mzz0+16.5*mww0-k3*mzz0      );
N4=((n+1)^q4-n^q4)*(      -10*mww0-k4*mww0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*myy(j)-k1*mxx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mzz(j)+16.5*mww(j)-k3*mzz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mww(j)-k4*mww(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*myy(j)-k1*mxx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzz(j)+16.5*mww(j)-k3*mzz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mww(j)-k4*mww(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*myy(j)-k1*mxx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mzz(j)+16.5*mww(j-td/h)-k3*mzz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mww(j-td/h)-k4*mww(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myy(j)-k1*mxx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzz(j)-2/19*myy(j)-1.08/19*sin(mxx(j))+0.061/19*sin(2*mxx(j))-k2*myy(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mzz(j)+16.5*mww(j-td/h)-k3*mzz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mww(j-td/h)-k4*mww(j)          );
end
mxx1(n+1)=mxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myy1(n+1)=myy0+h^q2*N2/(gamma(q2)*q2);
mzz1(n+1)=mzz0+h^q3*N3/(gamma(q3)*q3);
mww1(n+1)=mww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxx(n+2)=mxx0+h^q1*(       314*myy(n+1)-k1*mxx(n+1)      +M1)/gamma(q1+2);
myy(n+2)=myy0+h^q2*(       1/19*mzz(n+1)-2/19*myy(n+1)-1.08/19*sin(mxx(n+1))+0.061/19*sin(2*mxx(n+1))-k2*myy(n+1)       +M2)/gamma(q2+2);
mzz(n+2)=mzz0+h^q3*(       -2.5*mzz(n+1)+16.5*mww(n+1-td/h)-k3*mzz(n+1)    +M3)/gamma(q3+2);
mww(n+2)=mww0+h^q4*(      -10*mww(n+1-td/h)-k4*mww(n+1)    +M4)/gamma(q4+2);
end


q1=0.92;q2=0.92;q3=0.92;q4=0.92;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mxxx0=0.01;myyy0=0.01;mzzz0=0.01;mwww0=0.01;   %initial value
td=0.1;

mxxx(N+1)=[0];myyy(N+1)=[0];mzzz(N+1)=[0];mwww(N+1)=[0]; %efficiency need improve
mxxx1(N+1)=[0];myyy1(N+1)=[0];mzzz1(N+1)=[0];mwww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mxxx1(1)=mxxx0+h^q1*(   314*myyy0   -k1*mxxx0  )/(gamma(q1)*q1);
myyy1(1)=myyy0+h^q2*(   1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)   -k2*myyy0   )/(gamma(q2)*q2);
mzzz1(1)=mzzz0+h^q3*(   -2.5*mzzz0+16.5*mwww0   -k3*mzzz0         )/(gamma(q3)*q3);
mwww1(1)=mwww0+h^q4*(   -10*mwww0   -k4*mwww0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mxxx(1)=mxxx0+h^q1*((  314*myyy0-k1*mxxx0  )+q1*( 314*myyy0-k1*mxxx0  ))/gamma(q1+2);
myyy(1)=myyy0+h^q2*((  1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0 )+q2*(  1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0  ))/gamma(q2+2);
mzzz(1)=mzzz0+h^q3*((  -2.5*mzzz0+16.5*mwww0-k3*mzzz0  )+q3*(   -2.5*mzzz0+16.5*mwww0-k3*mzzz0   ))/gamma(q3+2);
mwww(1)=mwww0+h^q4*(( -10*mwww0-k4*mwww0  )+q4*(  -10*mwww0-k4*mwww0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyy0-k1*mxxx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mzzz0+16.5*mwww0-k3*mzzz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mwww0-k4*mwww0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*myyy0-k1*mxxx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mzzz0+16.5*mwww0-k3*mzzz0              );
N4=((n+1)^q4-n^q4)*(      -10*mwww0-k4*mwww0                     );
mxxx1(n+1)=mxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyy1(n+1)=myyy0+h^q2*N2/(gamma(q2)*q2);
mzzz1(n+1)=mzzz0+h^q3*N3/(gamma(q3)*q3);
mwww1(n+1)=mwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mxxx(n+2)=mxxx0+h^q1*(    314*myyy(n+1)-k1*mxxx(n+1)      +M1)/gamma(q1+2);
myyy(n+2)=myyy0+h^q2*(     1/19*mzzz(n+1)-2/19*myyy(n+1)-1.08/19*sin(mxxx(n+1))+0.061/19*sin(2*mxxx(n+1))-k2*myyy(n+1)      +M2)/gamma(q2+2);
mzzz(n+2)=mzzz0+h^q3*(    -2.5*mzzz(n+1)+16.5*mwww(n+1)-k3*mzzz(n+1)         +M3)/gamma(q3+2);
mwww(n+2)=mwww0+h^q4*(    -10*mwww(n+1)-k4*mwww(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyy0-k1*mxxx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mzzz0+16.5*mwww0-k3*mzzz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mwww0-k4*mwww0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*myyy0-k1*mxxx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mzzz0+16.5*mwww0-k3*mzzz0   );
N4=((n+1)^q4-n^q4)*(     -10*mwww0-k4*mwww0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*myyy(j)-k1*mxxx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mzzz(j)+16.5*mwww(j)-k3*mzzz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mwww(j)-k4*mwww(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyy(j)-k1*mxxx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzz(j)+16.5*mwww(j)-k3*mzzz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mwww(j)-k4*mwww(j)          );
end   
mxxx1(n+1)=mxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyy1(n+1)=myyy0+h^q2*N2/(gamma(q2)*q2);
mzzz1(n+1)=mzzz0+h^q3*N3/(gamma(q3)*q3);
mwww1(n+1)=mwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxx(n+2)=mxxx0+h^q1*(       314*myyy(n+1)-k1*mxxx(n+1)     +M1)/gamma(q1+2);
myyy(n+2)=myyy0+h^q2*(       1/19*mzzz(n+1)-2/19*myyy(n+1)-1.08/19*sin(mxxx(n+1))+0.061/19*sin(2*mxxx(n+1))-k2*myyy(n+1)        +M2)/gamma(q2+2);
mzzz(n+2)=mzzz0+h^q3*(       -2.5*mzzz(n+1)+16.5*mwww(n+1)-k3*mzzz(n+1)      +M3)/gamma(q3+2);
mwww(n+2)=mwww0+h^q4*(     -10*mwww(n+1)-k4*mwww(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*myyy0-k1*mxxx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mzzz0+16.5*mwww0-k3*mzzz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mwww0-k4*mwww0       );
N1=((n+1)^q1-n^q1)*(        314*myyy0-k1*mxxx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mzzz0-2/19*myyy0-1.08/19*sin(mxxx0)+0.061/19*sin(2*mxxx0)-k2*myyy0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mzzz0+16.5*mwww0-k3*mzzz0      );
N4=((n+1)^q4-n^q4)*(      -10*mwww0-k4*mwww0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*myyy(j)-k1*mxxx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mzzz(j)+16.5*mwww(j)-k3*mzzz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mwww(j)-k4*mwww(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*myyy(j)-k1*mxxx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzz(j)+16.5*mwww(j)-k3*mzzz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwww(j)-k4*mwww(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*myyy(j)-k1*mxxx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mzzz(j)+16.5*mwww(j-td/h)-k3*mzzz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mwww(j-td/h)-k4*mwww(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyy(j)-k1*mxxx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzz(j)-2/19*myyy(j)-1.08/19*sin(mxxx(j))+0.061/19*sin(2*mxxx(j))-k2*myyy(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mzzz(j)+16.5*mwww(j-td/h)-k3*mzzz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwww(j-td/h)-k4*mwww(j)          );
end
mxxx1(n+1)=mxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyy1(n+1)=myyy0+h^q2*N2/(gamma(q2)*q2);
mzzz1(n+1)=mzzz0+h^q3*N3/(gamma(q3)*q3);
mwww1(n+1)=mwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxx(n+2)=mxxx0+h^q1*(       314*myyy(n+1)-k1*mxxx(n+1)      +M1)/gamma(q1+2);
myyy(n+2)=myyy0+h^q2*(       1/19*mzzz(n+1)-2/19*myyy(n+1)-1.08/19*sin(mxxx(n+1))+0.061/19*sin(2*mxxx(n+1))-k2*myyy(n+1)       +M2)/gamma(q2+2);
mzzz(n+2)=mzzz0+h^q3*(       -2.5*mzzz(n+1)+16.5*mwww(n+1-td/h)-k3*mzzz(n+1)    +M3)/gamma(q3+2);
mwww(n+2)=mwww0+h^q4*(      -10*mwww(n+1-td/h)-k4*mwww(n+1)    +M4)/gamma(q4+2);
end


q1=0.94;q2=0.94;q3=0.94;q4=0.94;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mxxxx0=0.01;myyyy0=0.01;mzzzz0=0.01;mwwww0=0.01;   %initial value
td=0.1;

mxxxx(N+1)=[0];myyyy(N+1)=[0];mzzzz(N+1)=[0];mwwww(N+1)=[0]; %efficiency need improve
mxxxx1(N+1)=[0];myyyy1(N+1)=[0];mzzzz1(N+1)=[0];mwwww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%    replacement start  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 delay term replace)
mxxxx1(1)=mxxxx0+h^q1*(   314*myyyy0   -k1*mxxxx0  )/(gamma(q1)*q1);
myyyy1(1)=myyyy0+h^q2*(   1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)   -k2*myyyy0   )/(gamma(q2)*q2);
mzzzz1(1)=mzzzz0+h^q3*(   -2.5*mzzzz0+16.5*mwwww0   -k3*mzzzz0         )/(gamma(q3)*q3);
mwwww1(1)=mwwww0+h^q4*(   -10*mwwww0   -k4*mwwww0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) delay term replace)     f( x0 , y0 , z0 delay term replace)
mxxxx(1)=mxxxx0+h^q1*((  314*myyyy0-k1*mxxxx0  )+q1*( 314*myyyy0-k1*mxxxx0  ))/gamma(q1+2);
myyyy(1)=myyyy0+h^q2*((  1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0 )+q2*(  1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0  ))/gamma(q2+2);
mzzzz(1)=mzzzz0+h^q3*((  -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0  )+q3*(   -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0   ))/gamma(q3+2);
mwwww(1)=mwwww0+h^q4*(( -10*mwwww0-k4*mwwww0  )+q4*(  -10*mwwww0-k4*mwwww0  ))/gamma(q4+2);
%++++++++++++++++++++Calculate the first few delay terms have no value ++++++ make the delay terms direct=initial value +++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyyy0-k1*mxxxx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mwwww0-k4*mwwww0            );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(     314*myyyy0-k1*mxxxx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0              );
N4=((n+1)^q4-n^q4)*(      -10*mwwww0-k4*mwwww0                     );
mxxxx1(n+1)=mxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyy1(n+1)=myyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzz1(n+1)=mzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwww1(n+1)=mwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)
mxxxx(n+2)=mxxxx0+h^q1*(    314*myyyy(n+1)-k1*mxxxx(n+1)      +M1)/gamma(q1+2);
myyyy(n+2)=myyyy0+h^q2*(     1/19*mzzzz(n+1)-2/19*myyyy(n+1)-1.08/19*sin(mxxxx(n+1))+0.061/19*sin(2*mxxxx(n+1))-k2*myyyy(n+1)      +M2)/gamma(q2+2);
mzzzz(n+2)=mzzzz0+h^q3*(    -2.5*mzzzz(n+1)+16.5*mwwww(n+1)-k3*mzzzz(n+1)         +M3)/gamma(q3+2);
mwwww(n+2)=mwwww0+h^q4*(    -10*mwwww(n+1)-k4*mwwww(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyyy0-k1*mxxxx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mwwww0-k4*mwwww0           );
%                       f( x0 , y0 , z0 delay term replace )
N1=((n+1)^q1-n^q1)*(   314*myyyy0-k1*mxxxx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0   );
N4=((n+1)^q4-n^q4)*(     -10*mwwww0-k4*mwwww0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*myyyy(j)-k1*mxxxx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mzzzz(j)+16.5*mwwww(j)-k3*mzzzz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mwwww(j)-k4*mwwww(j)      );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyyy(j)-k1*mxxxx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzzz(j)+16.5*mwwww(j)-k3*mzzzz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mwwww(j)-k4*mwwww(j)          );
end   
mxxxx1(n+1)=mxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyy1(n+1)=myyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzz1(n+1)=mzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwww1(n+1)=mwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  delay term replace)

mxxxx(n+2)=mxxxx0+h^q1*(       314*myyyy(n+1)-k1*mxxxx(n+1)     +M1)/gamma(q1+2);
myyyy(n+2)=myyyy0+h^q2*(       1/19*mzzzz(n+1)-2/19*myyyy(n+1)-1.08/19*sin(mxxxx(n+1))+0.061/19*sin(2*mxxxx(n+1))-k2*myyyy(n+1)        +M2)/gamma(q2+2);
mzzzz(n+2)=mzzzz0+h^q3*(       -2.5*mzzzz(n+1)+16.5*mwwww(n+1)-k3*mzzzz(n+1)      +M3)/gamma(q3+2);
mwwww(n+2)=mwwww0+h^q4*(     -10*mwwww(n+1)-k4*mwwww(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%--------Calculate the case that the later delay term has a value ---- Some delay variables are replaced by the value of the preceding variable----------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 delay term replace)
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*myyyy0-k1*mxxxx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mwwww0-k4*mwwww0       );
N1=((n+1)^q1-n^q1)*(        314*myyyy0-k1*mxxxx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mzzzz0-2/19*myyyy0-1.08/19*sin(mxxxx0)+0.061/19*sin(2*mxxxx0)-k2*myyyy0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mzzzz0+16.5*mwwww0-k3*mzzzz0      );
N4=((n+1)^q4-n^q4)*(      -10*mwwww0-k4*mwwww0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*myyyy(j)-k1*mxxxx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mzzzz(j)+16.5*mwwww(j)-k3*mzzzz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mwwww(j)-k4*mwwww(j)     );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*myyyy(j)-k1*mxxxx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzzz(j)+16.5*mwwww(j)-k3*mzzzz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwwww(j)-k4*mwwww(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*myyyy(j)-k1*mxxxx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mzzzz(j)+16.5*mwwww(j-td/h)-k3*mzzzz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mwwww(j-td/h)-k4*mwwww(j)   );
      %    f(  x(j) , y(j) , z(j)  delay term replace)
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyyy(j)-k1*mxxxx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzzz(j)-2/19*myyyy(j)-1.08/19*sin(mxxxx(j))+0.061/19*sin(2*mxxxx(j))-k2*myyyy(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mzzzz(j)+16.5*mwwww(j-td/h)-k3*mzzzz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwwww(j-td/h)-k4*mwwww(j)          );
end
mxxxx1(n+1)=mxxxx0+h^q1*N1/(gamma(q1)*q1);                        %Estimation, Equation (22)
myyyy1(n+1)=myyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzz1(n+1)=mzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwww1(n+1)=mwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) 

mxxxx(n+2)=mxxxx0+h^q1*(       314*myyyy(n+1)-k1*mxxxx(n+1)      +M1)/gamma(q1+2);
myyyy(n+2)=myyyy0+h^q2*(       1/19*mzzzz(n+1)-2/19*myyyy(n+1)-1.08/19*sin(mxxxx(n+1))+0.061/19*sin(2*mxxxx(n+1))-k2*myyyy(n+1)       +M2)/gamma(q2+2);
mzzzz(n+2)=mzzzz0+h^q3*(       -2.5*mzzzz(n+1)+16.5*mwwww(n+1-td/h)-k3*mzzzz(n+1)    +M3)/gamma(q3+2);
mwwww(n+2)=mwwww0+h^q4*(      -10*mwwww(n+1-td/h)-k4*mwwww(n+1)    +M4)/gamma(q4+2);
end


q1=0.95;q2=0.95;q3=0.95;q4=0.95;
k1=45;k2=45;k3=6;k4=0;
h=0.02;N=999;           
mxxxxx0=0.01;myyyyy0=0.01;mzzzzz0=0.01;mwwwww0=0.01;   %
td=0.08;

mxxxxx(N+1)=[0];myyyyy(N+1)=[0];mzzzzz(N+1)=[0];mwwwww(N+1)=[0]; %
mxxxxx1(N+1)=[0];myyyyy1(N+1)=[0];mzzzzz1(N+1)=[0];mwwwww1(N+1)=[0];

%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%
%model： 
%  f( x0 , y0 , z0 
mxxxxx1(1)=mxxxxx0+h^q1*(   314*myyyyy0   -k1*mxxxxx0  )/(gamma(q1)*q1);
myyyyy1(1)=myyyyy0+h^q2*(   1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)   -k2*myyyyy0   )/(gamma(q2)*q2);
mzzzzz1(1)=mzzzzz0+h^q3*(   -2.5*mzzzzz0+16.5*mwwwww0   -k3*mzzzzz0         )/(gamma(q3)*q3);
mwwwww1(1)=mwwwww0+h^q4*(   -10*mwwwww0   -k4*mwwwww0         )/(gamma(q4)*q4);
%    f( x1(1) , y1(1) , z1(1) )     f( x0 , y0 , z0
mxxxxx(1)=mxxxxx0+h^q1*((  314*myyyyy0-k1*mxxxxx0  )+q1*( 314*myyyyy0-k1*mxxxxx0  ))/gamma(q1+2);
myyyyy(1)=myyyyy0+h^q2*((  1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0 )+q2*(  1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0  ))/gamma(q2+2);
mzzzzz(1)=mzzzzz0+h^q3*((  -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0  )+q3*(   -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0   ))/gamma(q3+2);
mwwwww(1)=mwwwww0+h^q4*(( -10*mwwwww0-k4*mwwwww0  )+q4*(  -10*mwwwww0-k4*mwwwww0  ))/gamma(q4+2);
%+++++++++++++++++++++++
n=0;
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyyyy0-k1*mxxxxx0             );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(   1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0         );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(   -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0          );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(   -10*mwwwww0-k4*mwwwww0            );
%                       f( x0 , y0 , z0
N1=((n+1)^q1-n^q1)*(     314*myyyyy0-k1*mxxxxx0              );
N2=((n+1)^q2-n^q2)*(       1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0                );
N3=((n+1)^q3-n^q3)*(       -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0              );
N4=((n+1)^q4-n^q4)*(      -10*mwwwww0-k4*mwwwww0                     );
mxxxxx1(n+1)=mxxxxx0+h^q1*N1/(gamma(q1)*q1);                        %22）
myyyyy1(n+1)=myyyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzzz1(n+1)=mzzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwwww1(n+1)=mwwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  
mxxxxx(n+2)=mxxxxx0+h^q1*(    314*myyyyy(n+1)-k1*mxxxxx(n+1)      +M1)/gamma(q1+2);
myyyyy(n+2)=myyyyy0+h^q2*(     1/19*mzzzzz(n+1)-2/19*myyyyy(n+1)-1.08/19*sin(mxxxxx(n+1))+0.061/19*sin(2*mxxxxx(n+1))-k2*myyyyy(n+1)      +M2)/gamma(q2+2);
mzzzzz(n+2)=mzzzzz0+h^q3*(    -2.5*mzzzzz(n+1)+16.5*mwwwww(n+1)-k3*mzzzzz(n+1)         +M3)/gamma(q3+2);
mwwwww(n+2)=mwwwww0+h^q4*(    -10*mwwwww(n+1)-k4*mwwwww(n+1)          +M4)/gamma(q4+2);
for n=1:td/h               %    f( x0 , y0 , z0
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(   314*myyyyy0-k1*mxxxxx0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0       );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0        );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(      -10*mwwwww0-k4*mwwwww0           );
%                       f( x0 , y0 , z0 )
N1=((n+1)^q1-n^q1)*(   314*myyyyy0-k1*mxxxxx0     );
N2=((n+1)^q2-n^q2)*(       1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0    );
N3=((n+1)^q3-n^q3)*(    -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0   );
N4=((n+1)^q4-n^q4)*(     -10*mwwwww0-k4*mwwwww0     );
for j=1:n
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    314*myyyyy(j)-k1*mxxxxx(j)         );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(   1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)    );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -2.5*mzzzzz(j)+16.5*mwwwww(j)-k3*mzzzzz(j)          );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(   -10*mwwwww(j)-k4*mwwwww(j)      );
      %    f(  x(j) , y(j) , z(j) )
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyyyy(j)-k1*mxxxxx(j)           );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)            );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzzzz(j)+16.5*mwwwww(j)-k3*mzzzzz(j)                );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(       -10*mwwwww(j)-k4*mwwwww(j)          );
end   
mxxxxx1(n+1)=mxxxxx0+h^q1*N1/(gamma(q1)*q1);                        %）
myyyyy1(n+1)=myyyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzzz1(n+1)=mzzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwwww1(n+1)=mwwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  

mxxxxx(n+2)=mxxxxx0+h^q1*(       314*myyyyy(n+1)-k1*mxxxxx(n+1)     +M1)/gamma(q1+2);
myyyyy(n+2)=myyyyy0+h^q2*(       1/19*mzzzzz(n+1)-2/19*myyyyy(n+1)-1.08/19*sin(mxxxxx(n+1))+0.061/19*sin(2*mxxxxx(n+1))-k2*myyyyy(n+1)        +M2)/gamma(q2+2);
mzzzzz(n+2)=mzzzzz0+h^q3*(       -2.5*mzzzzz(n+1)+16.5*mwwwww(n+1)-k3*mzzzzz(n+1)      +M3)/gamma(q3+2);
mwwwww(n+2)=mwwwww0+h^q4*(     -10*mwwwww(n+1)-k4*mwwwww(n+1)        +M4)/gamma(q4+2);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------

for n=td/h+1:N
    %                    f( x0 , y0 , z0 
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(      314*myyyyy0-k1*mxxxxx0         );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(        1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0          );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0     );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(       -10*mwwwww0-k4*mwwwww0       );
N1=((n+1)^q1-n^q1)*(        314*myyyyy0-k1*mxxxxx0              );
N2=((n+1)^q2-n^q2)*(        1/19*mzzzzz0-2/19*myyyyy0-1.08/19*sin(mxxxxx0)+0.061/19*sin(2*mxxxxx0)-k2*myyyyy0               );
N3=((n+1)^q3-n^q3)*(         -2.5*mzzzzz0+16.5*mwwwww0-k3*mzzzzz0      );
N4=((n+1)^q4-n^q4)*(      -10*mwwwww0-k4*mwwwww0          );
for j=1:td/h   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(       314*myyyyy(j)-k1*mxxxxx(j)              );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(       1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(        -2.5*mzzzzz(j)+16.5*mwwwww(j)-k3*mzzzzz(j)           );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(         -10*mwwwww(j)-k4*mwwwww(j)     );
      %    f(  x(j) , y(j) , z(j)  
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(      314*myyyyy(j)-k1*mxxxxx(j)              );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(       1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)               );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(       -2.5*mzzzzz(j)+16.5*mwwwww(j)-k3*mzzzzz(j)             );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwwwww(j)-k4*mwwwww(j)           );
end   
for j=td/h+1:n   
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     314*myyyyy(j)-k1*mxxxxx(j)          );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)          );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(       -2.5*mzzzzz(j)+16.5*mwwwww(j-td/h)-k3*mzzzzz(j)               );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(     -10*mwwwww(j-td/h)-k4*mwwwww(j)   );
      %    f(  x(j) , y(j) , z(j)  )
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(       314*myyyyy(j)-k1*mxxxxx(j)            );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(        1/19*mzzzzz(j)-2/19*myyyyy(j)-1.08/19*sin(mxxxxx(j))+0.061/19*sin(2*mxxxxx(j))-k2*myyyyy(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(        -2.5*mzzzzz(j)+16.5*mwwwww(j-td/h)-k3*mzzzzz(j)              );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(        -10*mwwwww(j-td/h)-k4*mwwwww(j)          );
end
mxxxxx1(n+1)=mxxxxx0+h^q1*N1/(gamma(q1)*q1);                        %（22）
myyyyy1(n+1)=myyyyy0+h^q2*N2/(gamma(q2)*q2);
mzzzzz1(n+1)=mzzzzz0+h^q3*N3/(gamma(q3)*q3);
mwwwww1(n+1)=mwwwww0+h^q4*N4/(gamma(q4)*q4);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  )

mxxxxx(n+2)=mxxxxx0+h^q1*(       314*myyyyy(n+1)-k1*mxxxxx(n+1)      +M1)/gamma(q1+2);
myyyyy(n+2)=myyyyy0+h^q2*(       1/19*mzzzzz(n+1)-2/19*myyyyy(n+1)-1.08/19*sin(mxxxxx(n+1))+0.061/19*sin(2*mxxxxx(n+1))-k2*myyyyy(n+1)       +M2)/gamma(q2+2);
mzzzzz(n+2)=mzzzzz0+h^q3*(       -2.5*mzzzzz(n+1)+16.5*mwwwww(n+1-td/h)-k3*mzzzzz(n+1)    +M3)/gamma(q3+2);
mwwwww(n+2)=mwwwww0+h^q4*(      -10*mwwwww(n+1-td/h)-k4*mwwwww(n+1)    +M4)/gamma(q4+2);
end
mde=8*m;mmde=8*mm;mmmde=8*mmm;mmmmde=8*mmmm;mmmmmde=8*mmmmm;pwde=1.5*pw;pwwde=1.5*pww;pwwwde=1.5*pwww;pwwwwde=1.5*pwwww;pwwwwwde=1.5*pwwwww;
mxde=0.7*mx;mxxde=0.7*mxx;mxxxde=0.7*mxxx;mxxxxde=0.7*mxxxx;mxxxxxde=0.7*mxxxxx;mzde=15*mz;mzzde=15*mzz;mzzzde=15*mzzz;mzzzzde=15*mzzzz;mzzzzzde=15*mzzzzz;
zo1=zeros(1,length(xx(1,:)))+0.9;
zo2=zeros(1,length(xx(1,:)))+0.8;
zo3=zeros(1,length(xx(1,:)))+0.7;
zo4=zeros(1,length(xx(1,:)))+0.6;
zo5=zeros(1,length(xx(1,:)))+0.5;

t=0:0.01:10;
figure(1)
plot3(zo1,t,w,'Linewidth',0.5);
hold on
plot3(zo1,t,px,'Linewidth',0.5);
hold on
plot3(zo1,t,mzde,'Linewidth',0.5);
hold on
plot3(zo2,t,ww,'Linewidth',0.5);
hold on
plot3(zo2,t,pxx,'Linewidth',0.5);
hold on
plot3(zo2,t,mzzde,'Linewidth',0.5);
hold on
plot3(zo3,t,www,'Linewidth',0.5);
hold on
plot3(zo3,t,pxxx,'Linewidth',0.5);
hold on
plot3(zo3,t,mzzzde,'Linewidth',0.5);
hold on
plot3(zo4,t,wwww,'Linewidth',0.5);
hold on
plot3(zo4,t,pxxxx,'Linewidth',0.5);
hold on
plot3(zo4,t,mzzzzde,'Linewidth',0.5);
hold on
plot3(zo5,t,wwwww,'Linewidth',0.5);
hold on
plot3(zo5,t,pxxxxx,'Linewidth',0.5);
hold on
plot3(zo5,t,mzzzzzde,'Linewidth',0.5);
xlabel('\alpha','fontsize',15);
ylabel('t','fontsize',15);
zlabel('x_4','fontsize',15);
legend('proposed method','PID control','existing method')
figure(2)
plot3(zo1,t,g,'Linewidth',0.5);
hold on
plot3(zo1,t,py,'Linewidth',0.5);
hold on
plot3(zo1,t,mw,'Linewidth',0.5);
hold on
plot3(zo2,t,gg,'Linewidth',0.5);
hold on
plot3(zo2,t,pyy,'Linewidth',0.5);
hold on
plot3(zo2,t,mww,'Linewidth',0.5);
hold on
plot3(zo3,t,ggg,'Linewidth',0.5);
hold on
plot3(zo3,t,pyyy,'Linewidth',0.5);
hold on
plot3(zo3,t,mwww,'Linewidth',0.5);
hold on
plot3(zo4,t,gggg,'Linewidth',0.5);
hold on
plot3(zo4,t,pyyyy,'Linewidth',0.5);
hold on
plot3(zo4,t,mwwww,'Linewidth',0.5);
hold on
plot3(zo5,t,ggggg,'Linewidth',0.5);
hold on
plot3(zo5,t,pyyyyy,'Linewidth',0.5);
hold on
plot3(zo5,t,mwwwww,'Linewidth',0.5);
xlabel('\alpha','fontsize',15);
ylabel('t','fontsize',15);
zlabel('x_5','fontsize',15);
legend('proposed method','PID control','existing method')
figure(3)
plot3(zo1,t,mde,'Linewidth',0.5);
hold on
plot3(zo1,t,pwde,'Linewidth',0.5);
hold on
plot3(zo1,t,mxde,'Linewidth',0.5);
hold on
plot3(zo2,t,mmde,'Linewidth',0.5);
hold on
plot3(zo2,t,pwwde,'Linewidth',0.5);
hold on
plot3(zo2,t,mxxde,'Linewidth',0.5);
hold on
plot3(zo3,t,mmmde,'Linewidth',0.5);
hold on
plot3(zo3,t,pwwwde,'Linewidth',0.5);
hold on
plot3(zo3,t,mxxxde,'Linewidth',0.5);
hold on
plot3(zo4,t,mmmmde,'Linewidth',0.5);
hold on
plot3(zo4,t,pwwwwde,'Linewidth',0.5);
hold on
plot3(zo4,t,mxxxxde,'Linewidth',0.5);
hold on
plot3(zo5,t,mmmmmde,'Linewidth',0.5);
hold on
plot3(zo5,t,pwwwwwde,'Linewidth',0.5);
hold on
plot3(zo5,t,mxxxxxde,'Linewidth',0.5);
xlabel('\alpha','fontsize',15);
ylabel('t','fontsize',15);
zlabel('x_6','fontsize',15);
legend('proposed method','PID control','existing method')
