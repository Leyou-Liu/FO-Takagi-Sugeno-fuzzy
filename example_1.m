clear;
clc;
%%              
A11=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 17231/16951 -2/9 0;0 0 0 0 0 0];
A12=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 1577/16951 -2/9 0;0 0 0 0 0 0];
Ad=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 -7/45;0 0 0 0 0 -10];
%B=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
%B=[0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];%%%             
B=eye(6);  %% 
I=eye(6);
setlmis([]);
X=lmivar(1,[6 1]);
M1=lmivar(2,[6 6]);  %Y1
M2=lmivar(2,[6 6]);  %Y2
Q=lmivar(1,[6 1]);

lmiterm([1 1 1 X],1,A11','s');
lmiterm([1 1 1 -M1],1,B','s');
lmiterm([1 1 2 0],Ad);
lmiterm([1 1 3 X],1,1);
lmiterm([1 2 2 Q],-1,1);
lmiterm([1 2 3 0],0);
lmiterm([1 3 3 Q],1,1);
lmiterm([1 3 3 0],-2*I);

lmiterm([2 1 1 X],1,A12','s');
lmiterm([2 1 1 -M2],1,B','s');
lmiterm([2 1 2 0],Ad);
lmiterm([2 1 3 X],1,1);
lmiterm([2 2 2 Q],-1,1);
lmiterm([2 2 3 0],0);
lmiterm([2 3 3 Q],1,1);
lmiterm([2 3 3 0],-2*I);

lmiterm([3 1 1 X],1,A11'+A12','s');
lmiterm([3 1 1 -M1],1,B','s');
lmiterm([3 1 1 -M2],1,B','s');
lmiterm([3 1 2 0],2*Ad);
lmiterm([3 1 3 X],2,1);
lmiterm([3 2 2 Q],-2,1);
lmiterm([3 2 3 0],0);
lmiterm([3 3 3 Q],2,1);
lmiterm([3 3 3 0],-4*I);

lmiterm([-4 1 1 X],1,1);
lmiterm([-5 1 1 Q],1,1);

lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
MATQ=dec2mat(lmisys,xfeas,X);
MATS=dec2mat(lmisys,xfeas,Q);
MATM1=dec2mat(lmisys,xfeas,M1);
MATM2=dec2mat(lmisys,xfeas,M2);
muk1=MATM1*inv(MATQ)
muk2=MATM2*inv(MATQ)
% 
% 
% Tya=0.1;
% A11=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 17231/16951 -2/9 0;0 0 0 0 0 0];
% A12=[0 1 0 0 0 0;0 0 1 0 0 0;-24 -24 -3 0 0 0;0 0 0 0 314 0;6.4 0 0.8 1577/16951 -2/9 0;0 0 0 0 0 0];
% Ad=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0;0 0 0 0 0 -7/45;0 0 0 0 0 -1/Tya];
% B=eye(6);
% % B=[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];%g
% % B=[0 0 0;0 0 0;0 0 0;1 0 0;0 1 0;0 0 1];
% I=eye(6);
% C=eye(6);
% % Bw=1.04*eye(6);
% Bw=1*eye(6);yip=0.08;ga=1;   %0.2 xt
% setlmis([]);
% X=lmivar(1,[6 1]);
% M1=lmivar(2,[6 6]);
% M2=lmivar(2,[6 6]);
% Q=lmivar(1,[6 1]);
% 
% lmiterm([1 1 1 X],1,A11','s');
% lmiterm([1 1 1 -M1],1,B','s');
% lmiterm([1 1 1 Q],1,1);
% lmiterm([1 1 1 0],C'*C);
% lmiterm([1 1 1 0],-I);
% 
% lmiterm([1 1 2 X],1,Ad);
% lmiterm([1 1 3 X],1,Bw);
% lmiterm([1 1 4 0],I);
% 
% lmiterm([1 2 2 Q],-1,1);
% lmiterm([1 2 3 0],0);
% lmiterm([1 2 4 0],0);
% lmiterm([1 3 3 0],-ga*ga*I);
% lmiterm([1 3 4 0],0);
% lmiterm([1 4 4 0],-I);
% 
% lmiterm([2 1 1 X],1,A12','s');
% lmiterm([2 1 1 -M2],1,B','s');
% lmiterm([2 1 1 Q],1,1);
% lmiterm([2 1 1 0],C'*C);
% lmiterm([2 1 1 0],-I);
% lmiterm([2 1 2 X],1,Ad);
% lmiterm([2 1 3 X],1,Bw);
% lmiterm([2 1 4 0],I);
% lmiterm([2 2 2 Q],-1,1);
% lmiterm([2 2 3 0],0);
% lmiterm([2 2 4 0],0);
% lmiterm([2 3 3 0],-ga*ga*I);
% lmiterm([2 3 4 0],0);
% lmiterm([2 4 4 0],-I);
% 
% lmiterm([3 1 1 X],1,A11'+A12','s');
% lmiterm([3 1 1 -M1],1,B','s');
% lmiterm([3 1 1 -M2],1,B','s');
% lmiterm([3 1 1 Q],1,1);
% lmiterm([3 1 1 0],C'*C);
% lmiterm([3 1 1 0],-I);
% lmiterm([3 1 2 X],1,Ad);
% lmiterm([3 1 3 X],1,Bw);
% lmiterm([3 1 4 0],I);
% lmiterm([3 2 2 Q],-1,1);
% lmiterm([3 2 3 0],0);
% lmiterm([3 2 4 0],0);
% lmiterm([3 3 3 0],-ga*ga*I);
% lmiterm([3 3 4 0],0);
% lmiterm([3 4 4 0],-I);
% 
% lmiterm([-4 1 1 X],1,1);
% lmiterm([-5 1 1 Q],1,1);
% 
% lmisys=getlmis;
% [tmin,xfeas]=feasp(lmisys);
% MATQ=dec2mat(lmisys,xfeas,X);
% MATS=dec2mat(lmisys,xfeas,Q);
% MATM1=dec2mat(lmisys,xfeas,M1);
% MATM2=dec2mat(lmisys,xfeas,M2);
% tmin;
% k1=MATM1/MATQ;
% k2=MATM2/MATQ;


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
%    -2.0080   -2.0222    1.9867    1.4769   -0.3993  -25.8605]; 

%
KKK1 =[0.0094    0.0034   -0.2119   -9.3270 -157.5028    0.0603
   -6.4301   -0.0310   -0.7939 -157.4898   -9.1058    0.0349
   -2.0080   -2.0222    1.9866    1.4774   -0.3522  -25.8593];
KKK2 =[0.0095    0.0034   -0.2136   -9.3270 -157.0453    0.0732
   -6.4297   -0.0309   -0.8043 -157.0236   -9.1047    0.1578
   -2.0080   -2.0222    1.9867    1.4769   -0.3993  -25.8605];


zhengtii=diag(0.6*(muk1)+0.4*(muk2));
zhengti=diag(zhengtii);
zhenggti=diag(zhenggti);

%%
q1=0.9;
q2=0.9;
q3=0.9;
q4=0.9;
q5=0.9;
q6=0.9;  


aw=0.0013;rw=0.0023;fp=0.052;Tya=0.1;   %   
aw=0.004;rw=0.00279;fp=0.052;Tya=0.1; %x    
% aw=0.0044;rw=0.0017;fp=0.052;Tya=0.1; %x     
h=0.01;N=500;           
v1=zhengti(1,1); v2=zhengti(1,2); v3=zhengti(1,3); v4=zhengti(1,4);v5=zhengti(1,5);v6=zhengti(1,6);
v7=zhengti(2,1); v8=zhengti(2,2); v9=zhengti(2,3); v10=zhengti(2,4);v11=zhengti(2,5);v12=zhengti(2,6);
v13=zhengti(3,1);v14=zhengti(3,2);v15=zhengti(3,3);v16=zhengti(3,4);v17=zhengti(3,5);v18=zhengti(3,6);
v19=zhengti(4,1);v20=zhengti(4,2);v21=zhengti(4,3);v22=zhengti(4,4);v23=zhengti(4,5);v24=zhengti(4,6);
v25=zhengti(5,1);v26=zhengti(5,2);v27=zhengti(5,3);v28=zhengti(5,4);v29=zhengti(5,5);v30=zhengti(5,6);
v31=zhengti(6,1);v32=zhengti(6,2);v33=zhengti(6,3);v34=zhengti(6,4);v35=zhengti(6,5);v36=zhengti(6,6);

xx0=0.01;
yy0=0.01;
zz0=0.01;
ww0=0.01;
gg0=0.01;
mm0=0.01;   %

td2=0.1;

xx(N+1)=[0];
yy(N+1)=[0];
zz(N+1)=[0];
ww(N+1)=[0];
gg(N+1)=[0];
mm(N+1)=[0];  %

xx1(N+1)=[0];
yy1(N+1)=[0];
zz1(N+1)=[0];
ww1(N+1)=[0];
gg1(N+1)=[0];
mm1(N+1)=[0];

%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
%： 
%  f( x0 , y0 , z0 )
xx1(1)=xx0+h^q1*(  yy0      +v1*xx0+v2*yy0+v3*zz0+v4*ww0+v5*gg0+v6*mm0   )/(gamma(q1)*q1);
yy1(1)=yy0+h^q2*(  zz0      +v7*xx0+v8*yy0+v9*zz0+v10*ww0+v11*gg0+v12*mm0   )/(gamma(q2)*q2);
zz1(1)=zz0+h^q3*(  -24*xx0-24*yy0-3*zz0+mm0        +v13*xx0+v14*yy0+v15*zz0+v16*ww0+v17*gg0+v18*mm0    )/(gamma(q3)*q3);
ww1(1)=ww0+h^q4*(  314*gg0    +v19*xx0+v20*yy0+v21*zz0+v22*ww0+v23*gg0+v24*mm0   )/(gamma(q4)*q4);
gg1(1)=gg0+h^q5*(  6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0   )/(gamma(q5)*q5);
mm1(1)=mm0+h^q6*(  -1/Tya*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0   )/(gamma(q6)*q6);
%    f( x1(1) , y1(1) , z1(1) )     f( x0 , y0 , z0 )
xx(1)=xx0+h^q1*((  yy0      +v1*xx0+v2*yy0+v3*zz0+v4*ww0+v5*gg0+v6*mm0  )+q1*(  yy0      +v1*xx0+v2*yy0+v3*zz0+v4*ww0+v5*gg0+v6*mm0    ))/gamma(q1+2);
yy(1)=yy0+h^q2*((    zz0      +v7*xx0+v8*yy0+v9*zz0+v10*ww0+v11*gg0+v12*mm0   )+q2*(  zz0      +v7*xx0+v8*yy0+v9*zz0+v10*ww0+v11*gg0+v12*mm0     ))/gamma(q2+2);
zz(1)=zz0+h^q3*((    -24*xx0-24*yy0-3*zz0+mm0        +v13*xx0+v14*yy0+v15*zz0+v16*ww0+v17*gg0+v18*mm0   )+q3*(  -24*xx0-24*yy0-3*zz0+mm0        +v13*xx0+v14*yy0+v15*zz0+v16*ww0+v17*gg0+v18*mm0   ))/gamma(q3+2);
ww(1)=ww0+h^q4*((   314*gg0    +v19*xx0+v20*yy0+v21*zz0+v22*ww0+v23*gg0+v24*mm0   )+q4*(  314*gg0    +v19*xx0+v20*yy0+v21*zz0+v22*ww0+v23*gg0+v24*mm0   ))/gamma(q4+2);
gg(1)=gg0+h^q5*((   6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0     )+q5*(  6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0   ))/gamma(q5+2);LL1';
mm(1)=mm0+h^q6*((   -1/Tya*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0  )+q6*(  -10*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0   ))/gamma(q6+2);
%++++++++++++++++++++++++
for n=0:td2/h-1               %    f( x0 , y0 , z0 )
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(     yy0      +v1*xx0+v2*yy0+v3*zz0+v4*ww0+v5*gg0+v6*mm0     );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     zz0      +v7*xx0+v8*yy0+v9*zz0+v10*ww0+v11*gg0+v12*mm0    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -24*xx0-24*yy0-3*zz0+mm0        +v13*xx0+v14*yy0+v15*zz0+v16*ww0+v17*gg0+v18*mm0  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*gg0    +v19*xx0+v20*yy0+v21*zz0+v22*ww0+v23*gg0+v24*mm0  );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(    6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0   );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(   -1/Tya*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0    );KKK1+KKK2;
%                       f( x0 , y0 , z0  )
N1=((n+1)^q1-n^q1)*(    yy0      +v1*xx0+v2*yy0+v3*zz0+v4*ww0+v5*gg0+v6*mm0  );
N2=((n+1)^q2-n^q2)*(     zz0      +v7*xx0+v8*yy0+v9*zz0+v10*ww0+v11*gg0+v12*mm0     );
N3=((n+1)^q3-n^q3)*(    -24*xx0-24*yy0-3*zz0+mm0        +v13*xx0+v14*yy0+v15*zz0+v16*ww0+v17*gg0+v18*mm0 );
N4=((n+1)^q4-n^q4)*(      314*gg0    +v19*xx0+v20*yy0+v21*zz0+v22*ww0+v23*gg0+v24*mm0     );
N5=((n+1)^q5-n^q5)*(   6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0   );
N6=((n+1)^q6-n^q6)*(     -1/Tya*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0       );
for j=1:n   %    f(  x(j) , y(j) , z(j)  )
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    yy(j)      +v1*xx(j)+v2*yy(j)+v3*zz(j)+v4*ww(j)+v5*gg(j)+v6*mm(j)   );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     zz(j)      +v7*xx(j)+v8*yy(j)+v9*zz(j)+v10*ww(j)+v11*gg(j)+v12*mm(j)     );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +v13*xx(j)+v14*yy(j)+v15*zz(j)+v16*ww(j)+v17*gg(j)+v18*mm(j)   );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aw*sin(fp*j)+rw*rand(1)+    314*gg(j)    +v19*xx(j)+v20*yy(j)+v21*zz(j)+v22*ww(j)+v23*gg(j)+v24*mm(j)   );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*( aw*sin(fp*j)+rw*rand(1)+   6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j)      );LL2';zhenggti';
        M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*( aw*sin(fp*j)+rw*rand(1)+  -1/Tya*mm(j)   +v31*xx(j)+v32*yy(j)+v33*zz(j)+v34*ww(j)+v35*gg(j)+v36*mm(j)    );
      %    f(  x(j) , y(j) , z(j)  )
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(      yy(j)      +v1*xx(j)+v2*yy(j)+v3*zz(j)+v4*ww(j)+v5*gg(j)+v6*mm(j)    );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(    zz(j)      +v7*xx(j)+v8*yy(j)+v9*zz(j)+v10*ww(j)+v11*gg(j)+v12*mm(j)  );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(    -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +v13*xx(j)+v14*yy(j)+v15*zz(j)+v16*ww(j)+v17*gg(j)+v18*mm(j)    );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(  aw*sin(fp*j)+rw*rand(1)+   314*gg(j)    +v19*xx(j)+v20*yy(j)+v21*zz(j)+v22*ww(j)+v23*gg(j)+v24*mm(j)    );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*( aw*sin(fp*j)+rw*rand(1)+   6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j)     );
        N6=N6+((n-j+1)^q6-(n-j)^q6)*( aw*sin(fp*j)+rw*rand(1)+   -1/Tya*mm(j)   +v31*xx(j)+v32*yy(j)+v33*zz(j)+v34*ww(j)+v35*gg(j)+v36*mm(j)      );
end   
xx1(n+1)=xx0+h^q1*N1/(gamma(q1)*q1);                        %）
yy1(n+1)=yy0+h^q2*N2/(gamma(q2)*q2);
zz1(n+1)=zz0+h^q3*N3/(gamma(q3)*q3);
ww1(n+1)=ww0+h^q4*N4/(gamma(q4)*q4);
gg1(n+1)=gg0+h^q5*N5/(gamma(q5)*q5);
mm1(n+1)=mm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  )
xx(n+2)=xx0+h^q1*(   yy(n+1)      +v1*xx(n+1)+v2*yy(n+1)+v3*zz(n+1)+v4*ww(n+1)+v5*gg(n+1)+v6*mm(n+1)   +M1)/gamma(q1+2);
yy(n+2)=yy0+h^q2*(    zz(n+1)      +v7*xx(n+1)+v8*yy(n+1)+v9*zz(n+1)+v10*ww(n+1)+v11*gg(n+1)+v12*mm(n+1)        +M2)/gamma(q2+2);
zz(n+2)=zz0+h^q3*(     -24*xx(n+1)-24*yy(n+1)-3*zz(n+1)+mm(n+1)        +v13*xx(n+1)+v14*yy(n+1)+v15*zz(n+1)+v16*ww(n+1)+v17*gg(n+1)+v18*mm(n+1)    +M3)/gamma(q3+2);
ww(n+2)=ww0+h^q4*(     314*gg(n+1)    +v19*xx(n+1)+v20*yy(n+1)+v21*zz(n+1)+v22*ww(n+1)+v23*gg(n+1)+v24*mm(n+1)        +M4)/gamma(q4+2);
gg(n+2)=gg0+h^q5*(     6.4*xx(n+1)+0.8*zz(n+1)-3/23*sin(ww(n+1))+320/30429*sin(2*ww(n+1))-2/9*gg(n+1)-7/45*mm(n+1)  +v25*xx(n+1)+v26*yy(n+1)+v27*zz(n+1)+v28*ww(n+1)+v29*gg(n+1)+v30*mm(n+1)    +M5)/gamma(q5+2);
mm(n+2)=mm0+h^q6*(     -1/Tya*mm(n+1)   +v31*xx(n+1)+v32*yy(n+1)+v33*zz(n+1)+v34*ww(n+1)+v35*gg(n+1)+v36*mm(n+1)      +M6)/gamma(q6+2);

end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%---------------
for n=td2/h:N
    %                    f( x0 , y0 , z0 )
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(    yy0      +v1*xx0+v2*yy0+v3*zz0+v4*ww0+v5*gg0+v6*mm0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(      zz0      +v7*xx0+v8*yy0+v9*zz0+v10*ww0+v11*gg0+v12*mm0     );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -24*xx0-24*yy0-3*zz0+mm0        +v13*xx0+v14*yy0+v15*zz0+v16*ww0+v17*gg0+v18*mm0   );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*gg0    +v19*xx0+v20*yy0+v21*zz0+v22*ww0+v23*gg0+v24*mm0   );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(     6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0  );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(      -1/Tya*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0    );

N1=((n+1)^q1-n^q1)*(     yy0      +v1*xx0+v2*yy0+v3*zz0+v4*ww0+v5*gg0+v6*mm0       );
N2=((n+1)^q2-n^q2)*(      zz0      +v7*xx0+v8*yy0+v9*zz0+v10*ww0+v11*gg0+v12*mm0      );
N3=((n+1)^q3-n^q3)*(     -24*xx0-24*yy0-3*zz0+mm0        +v13*xx0+v14*yy0+v15*zz0+v16*ww0+v17*gg0+v18*mm0  );
N4=((n+1)^q4-n^q4)*(    314*gg0    +v19*xx0+v20*yy0+v21*zz0+v22*ww0+v23*gg0+v24*mm0   );
N5=((n+1)^q5-n^q5)*(     6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0  );
N6=((n+1)^q6-n^q6)*(    -1/Tya*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0  );
for j=1:td2/h      
       %                                       f(  x(j) , y(j) , z(j) ) 
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(      yy(j)      +v1*xx(j)+v2*yy(j)+v3*zz(j)+v4*ww(j)+v5*gg(j)+v6*mm(j)      );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zz(j)      +v7*xx(j)+v8*yy(j)+v9*zz(j)+v10*ww(j)+v11*gg(j)+v12*mm(j)  );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +v13*xx(j)+v14*yy(j)+v15*zz(j)+v16*ww(j)+v17*gg(j)+v18*mm(j)  );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aw*sin(fp*j)+rw*rand(1)+     314*gg(j)    +v19*xx(j)+v20*yy(j)+v21*zz(j)+v22*ww(j)+v23*gg(j)+v24*mm(j)     );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(  aw*sin(fp*j)+rw*rand(1)+   6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j)  );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(  aw*sin(fp*j)+rw*rand(1)+     -1/Tya*mm(j)   +v31*xx(j)+v32*yy(j)+v33*zz(j)+v34*ww(j)+v35*gg(j)+v36*mm(j)     );
      %    f(  x(j) , y(j) , z(j)  )
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yy(j)      +v1*xx(j)+v2*yy(j)+v3*zz(j)+v4*ww(j)+v5*gg(j)+v6*mm(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zz(j)      +v7*xx(j)+v8*yy(j)+v9*zz(j)+v10*ww(j)+v11*gg(j)+v12*mm(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +v13*xx(j)+v14*yy(j)+v15*zz(j)+v16*ww(j)+v17*gg(j)+v18*mm(j)  );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(  aw*sin(fp*j)+rw*rand(1)+     314*gg(j)    +v19*xx(j)+v20*yy(j)+v21*zz(j)+v22*ww(j)+v23*gg(j)+v24*mm(j)       );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(  aw*sin(fp*j)+rw*rand(1)+     6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j) );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(  aw*sin(fp*j)+rw*rand(1)+      -1/Tya*mm(j)   +v31*xx(j)+v32*yy(j)+v33*zz(j)+v34*ww(j)+v35*gg(j)+v36*mm(j)       );
end   
for j=td2/h+1:n               
       %                    f(  x(j) , y(j) , z(j) )       
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     yy(j)      +v1*xx(j)+v2*yy(j)+v3*zz(j)+v4*ww(j)+v5*gg(j)+v6*mm(j)     );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zz(j)      +v7*xx(j)+v8*yy(j)+v9*zz(j)+v10*ww(j)+v11*gg(j)+v12*mm(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*xx(j)-24*yy(j)-3*zz(j)+mm(j-td2/h)        +v13*xx(j)+v14*yy(j)+v15*zz(j)+v16*ww(j)+v17*gg(j)+v18*mm(j) );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aw*sin(fp*j)+rw*rand(1)+    314*gg(j)    +v19*xx(j)+v20*yy(j)+v21*zz(j)+v22*ww(j)+v23*gg(j)+v24*mm(j)      );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(   aw*sin(fp*j)+rw*rand(1)+    6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j-td2/h)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j)   );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(    aw*sin(fp*j)+rw*rand(1)+   -1/Tya*mm(j-td2/h)   +v31*xx(j)+v32*yy(j)+v33*zz(j)+v34*ww(j)+v35*gg(j)+v36*mm(j)    );      
      %    f(  x(j) , y(j) , z(j)  )
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yy(j)      +v1*xx(j)+v2*yy(j)+v3*zz(j)+v4*ww(j)+v5*gg(j)+v6*mm(j)      );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zz(j)      +v7*xx(j)+v8*yy(j)+v9*zz(j)+v10*ww(j)+v11*gg(j)+v12*mm(j)   );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xx(j)-24*yy(j)-3*zz(j)+mm(j-td2/h)        +v13*xx(j)+v14*yy(j)+v15*zz(j)+v16*ww(j)+v17*gg(j)+v18*mm(j)   );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw*sin(fp*j)+rw*rand(1)+   314*gg(j)    +v19*xx(j)+v20*yy(j)+v21*zz(j)+v22*ww(j)+v23*gg(j)+v24*mm(j)   );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw*sin(fp*j)+rw*rand(1)+   6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j-td2/h)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j)  );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(  aw*sin(fp*j)+rw*rand(1)+     -1/Tya*mm(j-td2/h)   +v31*xx(j)+v32*yy(j)+v33*zz(j)+v34*ww(j)+v35*gg(j)+v36*mm(j)    );      
end
xx1(n+1)=xx0+h^q1*N1/(gamma(q1)*q1);                        %）
yy1(n+1)=yy0+h^q2*N2/(gamma(q2)*q2);
zz1(n+1)=zz0+h^q3*N3/(gamma(q2)*q3);
ww1(n+1)=ww0+h^q4*N4/(gamma(q2)*q4);
gg1(n+1)=gg0+h^q5*N5/(gamma(q5)*q5);
mm1(n+1)=mm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  )
xx(n+1)=xx0+h^q1*(     yy(n+1)      +v1*xx(n+1)+v2*yy(n+1)+v3*zz(n+1)+v4*ww(n+1)+v5*gg(n+1)+v6*mm(n+1)   +M1)/gamma(q1+2);
yy(n+1)=yy0+h^q2*(     zz(n+1)      +v7*xx(n+1)+v8*yy(n+1)+v9*zz(n+1)+v10*ww(n+1)+v11*gg(n+1)+v12*mm(n+1)     +M2)/gamma(q2+2);
zz(n+1)=zz0+h^q3*(     -24*xx(n+1)-24*yy(n+1)-3*zz(n+1)+mm(n+1-td2/h)        +v13*xx(n+1)+v14*yy(n+1)+v15*zz(n+1)+v16*ww(n+1)+v17*gg(n+1)+v18*mm(n+1)    +M3)/gamma(q3+2);
ww(n+1)=ww0+h^q4*(     314*gg(n+1)    +v19*xx(n+1)+v20*yy(n+1)+v21*zz(n+1)+v22*ww(n+1)+v23*gg(n+1)+v24*mm(n+1)      +M4)/gamma(q4+2);
gg(n+1)=gg0+h^q5*(     6.4*xx(n+1)+0.8*zz(n+1)-3/23*sin(ww(n+1))+320/30429*sin(2*ww(n+1))-2/9*gg(n+1)-7/45*mm(n+1-td2/h)  +v25*xx(n+1)+v26*yy(n+1)+v27*zz(n+1)+v28*ww(n+1)+v29*gg(n+1)+v30*mm(n+1)    +M5)/gamma(q5+2);
mm(n+1)=mm0+h^q6*(     -1/Tya*mm(n+1-td2/h)   +v31*xx(n+1)+v32*yy(n+1)+v33*zz(n+1)+v34*ww(n+1)+v35*gg(n+1)+v36*mm(n+1)      +M6)/gamma(q6+2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=0.01;y0=0.01;z0=0.01;w0=0.01;g0=0.01;m0=0.01;  
aaw=0;rrw=0;
x(N+1)=[0];y(N+1)=[0];z(N+1)=[0];w(N+1)=[0];g(N+1)=[0];m(N+1)=[0];  %
x1(N+1)=[0];y1(N+1)=[0];z1(N+1)=[0];w1(N+1)=[0];g1(N+1)=[0];m1(N+1)=[0];
%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%
%： 
%  f( x0 , y0 , z0)
x1(1)=x0+h^q1*(  y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0   )/(gamma(q1)*q1);
y1(1)=y0+h^q2*(  z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0   )/(gamma(q2)*q2);
z1(1)=z0+h^q3*(  -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0    )/(gamma(q3)*q3);
w1(1)=w0+h^q4*(  314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   )/(gamma(q4)*q4);
g1(1)=g0+h^q5*(  6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0   )/(gamma(q5)*q5);
m1(1)=m0+h^q6*(  -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0   )/(gamma(q6)*q6);
%    f( x1(1) , y1(1) , z1(1) )     f( x0 , y0 , z0 )
x(1)=x0+h^q1*((  y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0  )+q1*(  y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0    ))/gamma(q1+2);
y(1)=y0+h^q2*((    z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0   )+q2*(  z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0     ))/gamma(q2+2);
z(1)=z0+h^q3*((    -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0   )+q3*(  -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0   ))/gamma(q3+2);
w(1)=w0+h^q4*((   314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   )+q4*(  314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   ))/gamma(q4+2);
g(1)=g0+h^q5*((   6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0     )+q5*(  6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0   ))/gamma(q5+2);
m(1)=m0+h^q6*((   -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0  )+q6*(  -10*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0   ))/gamma(q6+2);
%++++++++++++++++++++++++
for n=0:td2/h-1               %    f( x0 , y0 , z0 )
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(     y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0     );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0  );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(    6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0   );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(   -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0    );
%                       f( x0 , y0 , z0 )
N1=((n+1)^q1-n^q1)*(    y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0  );
N2=((n+1)^q2-n^q2)*(     z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0     );
N3=((n+1)^q3-n^q3)*(    -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0 );
N4=((n+1)^q4-n^q4)*(      314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0     );
N5=((n+1)^q5-n^q5)*(   6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0   );
N6=((n+1)^q6-n^q6)*(     -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0       );
for j=1:n   %    f(  x(j) , y(j) , z(j)  )
M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(    y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)   );
 M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(     z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)     );
  M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*x(j)-24*y(j)-3*z(j)+m(j)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)   );
    M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*( aaw*sin(fp*j)+rrw*rand(1)+    314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)   );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*( aaw*sin(fp*j)+rrw*rand(1)+   6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)      );
        M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*( aaw*sin(fp*j)+rrw*rand(1)+  -1/Tya*m(j)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)    );
      %    f(  x(j) , y(j) , z(j)  )
 N1=N1+((n-j+1)^q1-(n-j)^q1)*(      y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)    );
  N2=N2+((n-j+1)^q2-(n-j)^q2)*(    z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)  );
   N3=N3+((n-j+1)^q3-(n-j)^q3)*(    -24*x(j)-24*y(j)-3*z(j)+m(j)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)    );
    N4=N4+((n-j+1)^q4-(n-j)^q4)*(  aaw*sin(fp*j)+rrw*rand(1)+   314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)    );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*( aaw*sin(fp*j)+rrw*rand(1)+   6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)     );
        N6=N6+((n-j+1)^q6-(n-j)^q6)*( aaw*sin(fp*j)+rrw*rand(1)+   -1/Tya*m(j)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)      );
end   
x1(n+1)=x0+h^q1*N1/(gamma(q1)*q1);                        %
y1(n+1)=y0+h^q2*N2/(gamma(q2)*q2);
z1(n+1)=z0+h^q3*N3/(gamma(q3)*q3);
w1(n+1)=w0+h^q4*N4/(gamma(q4)*q4);
g1(n+1)=g0+h^q5*N5/(gamma(q5)*q5);
m1(n+1)=m0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) )
x(n+2)=x0+h^q1*(   y(n+1)      +v1*x(n+1)+v2*y(n+1)+v3*z(n+1)+v4*w(n+1)+v5*g(n+1)+v6*m(n+1)   +M1)/gamma(q1+2);
y(n+2)=y0+h^q2*(    z(n+1)      +v7*x(n+1)+v8*y(n+1)+v9*z(n+1)+v10*w(n+1)+v11*g(n+1)+v12*m(n+1)        +M2)/gamma(q2+2);
z(n+2)=z0+h^q3*(     -24*x(n+1)-24*y(n+1)-3*z(n+1)+m(n+1)        +v13*x(n+1)+v14*y(n+1)+v15*z(n+1)+v16*w(n+1)+v17*g(n+1)+v18*m(n+1)    +M3)/gamma(q3+2);
w(n+2)=w0+h^q4*(     314*g(n+1)    +v19*x(n+1)+v20*y(n+1)+v21*z(n+1)+v22*w(n+1)+v23*g(n+1)+v24*m(n+1)        +M4)/gamma(q4+2);
g(n+2)=g0+h^q5*(     6.4*x(n+1)+0.8*z(n+1)-3/23*sin(w(n+1))+320/30429*sin(2*w(n+1))-2/9*g(n+1)-7/45*m(n+1)  +v25*x(n+1)+v26*y(n+1)+v27*z(n+1)+v28*w(n+1)+v29*g(n+1)+v30*m(n+1)    +M5)/gamma(q5+2);
m(n+2)=m0+h^q6*(     -1/Tya*m(n+1)   +v31*x(n+1)+v32*y(n+1)+v33*z(n+1)+v34*w(n+1)+v35*g(n+1)+v36*m(n+1)      +M6)/gamma(q6+2);

end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------
for n=td2/h:N
    %                    f( x0 , y0 , z0 )
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(    y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0    );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(      z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0     );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(     -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0   );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(     6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0  );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(      -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0    );

N1=((n+1)^q1-n^q1)*(     y0      +v1*x0+v2*y0+v3*z0+v4*w0+v5*g0+v6*m0       );
N2=((n+1)^q2-n^q2)*(      z0      +v7*x0+v8*y0+v9*z0+v10*w0+v11*g0+v12*m0      );
N3=((n+1)^q3-n^q3)*(     -24*x0-24*y0-3*z0+m0        +v13*x0+v14*y0+v15*z0+v16*w0+v17*g0+v18*m0  );
N4=((n+1)^q4-n^q4)*(    314*g0    +v19*x0+v20*y0+v21*z0+v22*w0+v23*g0+v24*m0   );
N5=((n+1)^q5-n^q5)*(     6.4*x0+0.8*z0-3/23*sin(w0)+320/30429*sin(2*w0)-2/9*g0-7/45*m0  +v25*x0+v26*y0+v27*z0+v28*w0+v29*g0+v30*m0  );
N6=((n+1)^q6-n^q6)*(    -1/Tya*m0   +v31*x0+v32*y0+v33*z0+v34*w0+v35*g0+v36*m0  );
for j=1:td2/h      
       %                                       f(  x(j) , y(j) , z(j)  ) 
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(      y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)      );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)  );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -24*x(j)-24*y(j)-3*z(j)+m(j)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)  );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aaw*sin(fp*j)+rrw*rand(1)+     314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)     );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(  aaw*sin(fp*j)+rrw*rand(1)+   6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)  );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(  aaw*sin(fp*j)+rrw*rand(1)+     -1/Tya*m(j)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)     );
      %    f(  x(j) , y(j) , z(j)  )
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)       );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)      );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*x(j)-24*y(j)-3*z(j)+m(j)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)  );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(  aaw*sin(fp*j)+rrw*rand(1)+     314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)       );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(  aaw*sin(fp*j)+rrw*rand(1)+     6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j) );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(  aaw*sin(fp*j)+rrw*rand(1)+      -1/Tya*m(j)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)       );
end   
for j=td2/h+1:n               
       %                    f(  x(j) , y(j) , z(j)  )       
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(     y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)     );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)   );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(      -24*x(j)-24*y(j)-3*z(j)+m(j-td2/h)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j) );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aaw*sin(fp*j)+rrw*rand(1)+    314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)      );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(   aaw*sin(fp*j)+rrw*rand(1)+    6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j-td2/h)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)   );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(    aaw*sin(fp*j)+rrw*rand(1)+   -1/Tya*m(j-td2/h)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)    );      
      %    f(  x(j) , y(j) , z(j)  )
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     y(j)      +v1*x(j)+v2*y(j)+v3*z(j)+v4*w(j)+v5*g(j)+v6*m(j)      );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      z(j)      +v7*x(j)+v8*y(j)+v9*z(j)+v10*w(j)+v11*g(j)+v12*m(j)   );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*x(j)-24*y(j)-3*z(j)+m(j-td2/h)        +v13*x(j)+v14*y(j)+v15*z(j)+v16*w(j)+v17*g(j)+v18*m(j)   );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aaw*sin(fp*j)+rrw*rand(1)+   314*g(j)    +v19*x(j)+v20*y(j)+v21*z(j)+v22*w(j)+v23*g(j)+v24*m(j)   );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aaw*sin(fp*j)+rrw*rand(1)+   6.4*x(j)+0.8*z(j)-3/23*sin(w(j))+320/30429*sin(2*w(j))-2/9*g(j)-7/45*m(j-td2/h)  +v25*x(j)+v26*y(j)+v27*z(j)+v28*w(j)+v29*g(j)+v30*m(j)  );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(  aaw*sin(fp*j)+rrw*rand(1)+     -1/Tya*m(j-td2/h)   +v31*x(j)+v32*y(j)+v33*z(j)+v34*w(j)+v35*g(j)+v36*m(j)    );      
end
x1(n+1)=x0+h^q1*N1/(gamma(q1)*q1);                        %）
y1(n+1)=y0+h^q2*N2/(gamma(q2)*q2);
z1(n+1)=z0+h^q3*N3/(gamma(q2)*q3);
w1(n+1)=w0+h^q4*N4/(gamma(q2)*q4);
g1(n+1)=g0+h^q5*N5/(gamma(q5)*q5);
m1(n+1)=m0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1)  )
x(n+1)=x0+h^q1*(     y(n+1)      +v1*x(n+1)+v2*y(n+1)+v3*z(n+1)+v4*w(n+1)+v5*g(n+1)+v6*m(n+1)   +M1)/gamma(q1+2);
y(n+1)=y0+h^q2*(     z(n+1)      +v7*x(n+1)+v8*y(n+1)+v9*z(n+1)+v10*w(n+1)+v11*g(n+1)+v12*m(n+1)     +M2)/gamma(q2+2);
z(n+1)=z0+h^q3*(     -24*x(n+1)-24*y(n+1)-3*z(n+1)+m(n+1-td2/h)        +v13*x(n+1)+v14*y(n+1)+v15*z(n+1)+v16*w(n+1)+v17*g(n+1)+v18*m(n+1)    +M3)/gamma(q3+2);
w(n+1)=w0+h^q4*(     314*g(n+1)    +v19*x(n+1)+v20*y(n+1)+v21*z(n+1)+v22*w(n+1)+v23*g(n+1)+v24*m(n+1)      +M4)/gamma(q4+2);
g(n+1)=g0+h^q5*(     6.4*x(n+1)+0.8*z(n+1)-3/23*sin(w(n+1))+320/30429*sin(2*w(n+1))-2/9*g(n+1)-7/45*m(n+1-td2/h)  +v25*x(n+1)+v26*y(n+1)+v27*z(n+1)+v28*w(n+1)+v29*g(n+1)+v30*m(n+1)    +M5)/gamma(q5+2);
m(n+1)=m0+h^q6*(     -1/Tya*m(n+1-td2/h)   +v31*x(n+1)+v32*y(n+1)+v33*z(n+1)+v34*w(n+1)+v35*g(n+1)+v36*m(n+1)      +M6)/gamma(q6+2);
end
t=0:0.01:5;
figure(1) 
plot(t,ww*1.0,'Linewidth',0.4);
hold on 
plot(t,w*1,'Linewidth',0.4);  
xlabel('t/s','fontsize',15);
ylabel('{x_4}','fontsize',15);
legend('controller without u_z(t)','controller with u_z(t)')
figure(2)
plot(t,gg*1.,'Linewidth',0.4);
plot(t,gg*1.1,'Linewidth',0.4);%zhen
hold on
plot(t,g*1.0,'Linewidth',0.4);
plot(t,g*1.,'Linewidth',0.4);
xlabel('t/s','fontsize',15);
ylabel('{x_5}','fontsize',15);
legend('controller without u_z(t)','controller with u_z(t)')%10
figure(3)
plot(t,mm*1,'Linewidth',0.4);  
hold on
plot(t,m,'Linewidth',0.4); 
xlabel('t/s','fontsize',15);
ylabel('{x_6}','fontsize',15);
legend('controller without u_z(t)','controller with u_z(t)')
