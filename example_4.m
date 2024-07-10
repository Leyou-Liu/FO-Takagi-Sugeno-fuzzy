function example_4
zhengti=zeros(6);
v1=zhengti(1,1); v2=zhengti(1,2); v3=zhengti(1,3); v4=zhengti(1,4);v5=zhengti(1,5);v6=zhengti(1,6);
v7=zhengti(2,1); v8=zhengti(2,2); v9=zhengti(2,3); v10=zhengti(2,4);v11=zhengti(2,5);v12=zhengti(2,6);
v13=zhengti(3,1);v14=zhengti(3,2);v15=zhengti(3,3);v16=zhengti(3,4);v17=zhengti(3,5);v18=zhengti(3,6);
v19=zhengti(4,1);v20=zhengti(4,2);v21=zhengti(4,3);v22=zhengti(4,4);v23=zhengti(4,5);v24=zhengti(4,6);
v25=zhengti(5,1);v26=zhengti(5,2);v27=zhengti(5,3);v28=zhengti(5,4);v29=zhengti(5,5);v30=zhengti(5,6);
v31=zhengti(6,1);v32=zhengti(6,2);v33=zhengti(6,3);v34=zhengti(6,4);v35=zhengti(6,5);v36=zhengti(6,6);
q1=0.9;q2=0.9;q3=0.9;q4=0.9;q5=0.9;q6=0.9;  
% q1=0.98;q2=0.98;q3=0.98;q4=0.98;q5=0.98;q6=0.98;  
% q1=0.8;q2=0.8;q3=0.8;q4=0.8;q5=0.8;q6=0.8; 
aw=0.02;rw=0.012;fp=5; Tya=0.1; td2=0.1;  %0.1 0.02 0.012
h=0.01;N=2000;           
xx0=0.01;yy0=0.01;zz0=0.01;ww0=0.01;gg0=0.01;mm0=0.01;   %
xx(N+1)=[0];yy(N+1)=[0];zz(N+1)=[0];ww(N+1)=[0];gg(N+1)=[0];mm(N+1)=[0];  %
xx1(N+1)=[0];yy1(N+1)=[0];zz1(N+1)=[0];ww1(N+1)=[0];gg1(N+1)=[0];mm1(N+1)=[0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
gg(1)=gg0+h^q5*((   6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0     )+q5*(  6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0   ))/gamma(q5+2);
mm(1)=mm0+h^q6*((   -1/Tya*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0  )+q6*(  -10*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0   ))/gamma(q6+2);
%++++++++++++++++++++++++++
for n=0:td2/h-1               %    f( x0 , y0 , z0 
M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(     yy0      +v1*xx0+v2*yy0+v3*zz0+v4*ww0+v5*gg0+v6*mm0     );
M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(     zz0      +v7*xx0+v8*yy0+v9*zz0+v10*ww0+v11*gg0+v12*mm0    );
M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(      -24*xx0-24*yy0-3*zz0+mm0        +v13*xx0+v14*yy0+v15*zz0+v16*ww0+v17*gg0+v18*mm0  );
M4=(n^(q4+1)-(n-q4)*(n+1)^q4)*(     314*gg0    +v19*xx0+v20*yy0+v21*zz0+v22*ww0+v23*gg0+v24*mm0  );
M5=(n^(q5+1)-(n-q5)*(n+1)^q5)*(    6.4*xx0+0.8*zz0-3/23*sin(ww0)+320/30429*sin(2*ww0)-2/9*gg0-7/45*mm0  +v25*xx0+v26*yy0+v27*zz0+v28*ww0+v29*gg0+v30*mm0   );
M6=(n^(q6+1)-(n-q6)*(n+1)^q6)*(   -1/Tya*mm0   +v31*xx0+v32*yy0+v33*zz0+v34*ww0+v35*gg0+v36*mm0    );
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
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*( aw*sin(fp*j)+rw*rand(1)+   6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j)      );
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
%    f(  x1(n+1) , y1(n+1) , z1(n+1) 
xx(n+2)=xx0+h^q1*(   yy(n+1)      +v1*xx(n+1)+v2*yy(n+1)+v3*zz(n+1)+v4*ww(n+1)+v5*gg(n+1)+v6*mm(n+1)   +M1)/gamma(q1+2);
yy(n+2)=yy0+h^q2*(    zz(n+1)      +v7*xx(n+1)+v8*yy(n+1)+v9*zz(n+1)+v10*ww(n+1)+v11*gg(n+1)+v12*mm(n+1)        +M2)/gamma(q2+2);
zz(n+2)=zz0+h^q3*(     -24*xx(n+1)-24*yy(n+1)-3*zz(n+1)+mm(n+1)        +v13*xx(n+1)+v14*yy(n+1)+v15*zz(n+1)+v16*ww(n+1)+v17*gg(n+1)+v18*mm(n+1)    +M3)/gamma(q3+2);
ww(n+2)=ww0+h^q4*(     314*gg(n+1)    +v19*xx(n+1)+v20*yy(n+1)+v21*zz(n+1)+v22*ww(n+1)+v23*gg(n+1)+v24*mm(n+1)        +M4)/gamma(q4+2);
gg(n+2)=gg0+h^q5*(     6.4*xx(n+1)+0.8*zz(n+1)-3/23*sin(ww(n+1))+320/30429*sin(2*ww(n+1))-2/9*gg(n+1)-7/45*mm(n+1)  +v25*xx(n+1)+v26*yy(n+1)+v27*zz(n+1)+v28*ww(n+1)+v29*gg(n+1)+v30*mm(n+1)    +M5)/gamma(q5+2);
mm(n+2)=mm0+h^q6*(     -1/Tya*mm(n+1)   +v31*xx(n+1)+v32*yy(n+1)+v33*zz(n+1)+v34*ww(n+1)+v35*gg(n+1)+v36*mm(n+1)      +M6)/gamma(q6+2);

end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%----------------
for n=td2/h:N
    %                    f( x0 , y0 , z0
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
       %                                       f(  x(j) , y(j) , z(j)  ) 
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(      yy(j)      +v1*xx(j)+v2*yy(j)+v3*zz(j)+v4*ww(j)+v5*gg(j)+v6*mm(j)      );
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(      zz(j)      +v7*xx(j)+v8*yy(j)+v9*zz(j)+v10*ww(j)+v11*gg(j)+v12*mm(j)  );
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(     -24*xx(j)-24*yy(j)-3*zz(j)+mm(j)        +v13*xx(j)+v14*yy(j)+v15*zz(j)+v16*ww(j)+v17*gg(j)+v18*mm(j)  );
      M4=M4+((n-j+2)^(q4+1)+(n-j)^(q4+1)-2*(n-j+1)^(q4+1))*(  aw*sin(fp*j)+rw*rand(1)+     314*gg(j)    +v19*xx(j)+v20*yy(j)+v21*zz(j)+v22*ww(j)+v23*gg(j)+v24*mm(j)     );
      M5=M5+((n-j+2)^(q5+1)+(n-j)^(q5+1)-2*(n-j+1)^(q5+1))*(  aw*sin(fp*j)+rw*rand(1)+   6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j)  );
      M6=M6+((n-j+2)^(q6+1)+(n-j)^(q6+1)-2*(n-j+1)^(q6+1))*(  aw*sin(fp*j)+rw*rand(1)+     -1/Tya*mm(j)   +v31*xx(j)+v32*yy(j)+v33*zz(j)+v34*ww(j)+v35*gg(j)+v36*mm(j)     );
      %    f(  x(j) , y(j) , z(j)  
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
      %    f(  x(j) , y(j) , z(j)  
      N1=N1+((n-j+1)^q1-(n-j)^q1)*(     yy(j)      +v1*xx(j)+v2*yy(j)+v3*zz(j)+v4*ww(j)+v5*gg(j)+v6*mm(j)      );
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(      zz(j)      +v7*xx(j)+v8*yy(j)+v9*zz(j)+v10*ww(j)+v11*gg(j)+v12*mm(j)   );
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(     -24*xx(j)-24*yy(j)-3*zz(j)+mm(j-td2/h)        +v13*xx(j)+v14*yy(j)+v15*zz(j)+v16*ww(j)+v17*gg(j)+v18*mm(j)   );
      N4=N4+((n-j+1)^q4-(n-j)^q4)*(   aw*sin(fp*j)+rw*rand(1)+   314*gg(j)    +v19*xx(j)+v20*yy(j)+v21*zz(j)+v22*ww(j)+v23*gg(j)+v24*mm(j)   );
      N5=N5+((n-j+1)^q5-(n-j)^q5)*(   aw*sin(fp*j)+rw*rand(1)+   6.4*xx(j)+0.8*zz(j)-3/23*sin(ww(j))+320/30429*sin(2*ww(j))-2/9*gg(j)-7/45*mm(j-td2/h)  +v25*xx(j)+v26*yy(j)+v27*zz(j)+v28*ww(j)+v29*gg(j)+v30*mm(j)  );
      N6=N6+((n-j+1)^q6-(n-j)^q6)*(  aw*sin(fp*j)+rw*rand(1)+     -1/Tya*mm(j-td2/h)   +v31*xx(j)+v32*yy(j)+v33*zz(j)+v34*ww(j)+v35*gg(j)+v36*mm(j)    );      
end
xx1(n+1)=xx0+h^q1*N1/(gamma(q1)*q1);                        %2）
yy1(n+1)=yy0+h^q2*N2/(gamma(q2)*q2);
zz1(n+1)=zz0+h^q3*N3/(gamma(q2)*q3);
ww1(n+1)=ww0+h^q4*N4/(gamma(q2)*q4);
gg1(n+1)=gg0+h^q5*N5/(gamma(q5)*q5);
mm1(n+1)=mm0+h^q6*N6/(gamma(q6)*q6);
%    f(  x1(n+1) , y1(n+1) , z1(n+1) 
xx(n+1)=xx0+h^q1*(     yy(n+1)      +v1*xx(n+1)+v2*yy(n+1)+v3*zz(n+1)+v4*ww(n+1)+v5*gg(n+1)+v6*mm(n+1)   +M1)/gamma(q1+2);
yy(n+1)=yy0+h^q2*(     zz(n+1)      +v7*xx(n+1)+v8*yy(n+1)+v9*zz(n+1)+v10*ww(n+1)+v11*gg(n+1)+v12*mm(n+1)     +M2)/gamma(q2+2);
zz(n+1)=zz0+h^q3*(     -24*xx(n+1)-24*yy(n+1)-3*zz(n+1)+mm(n+1-td2/h)        +v13*xx(n+1)+v14*yy(n+1)+v15*zz(n+1)+v16*ww(n+1)+v17*gg(n+1)+v18*mm(n+1)    +M3)/gamma(q3+2);
ww(n+1)=ww0+h^q4*(     314*gg(n+1)    +v19*xx(n+1)+v20*yy(n+1)+v21*zz(n+1)+v22*ww(n+1)+v23*gg(n+1)+v24*mm(n+1)      +M4)/gamma(q4+2);
gg(n+1)=gg0+h^q5*(     6.4*xx(n+1)+0.8*zz(n+1)-3/23*sin(ww(n+1))+320/30429*sin(2*ww(n+1))-2/9*gg(n+1)-7/45*mm(n+1-td2/h)  +v25*xx(n+1)+v26*yy(n+1)+v27*zz(n+1)+v28*ww(n+1)+v29*gg(n+1)+v30*mm(n+1)    +M5)/gamma(q5+2);
mm(n+1)=mm0+h^q6*(     -1/Tya*mm(n+1-td2/h)   +v31*xx(n+1)+v32*yy(n+1)+v33*zz(n+1)+v34*ww(n+1)+v35*gg(n+1)+v36*mm(n+1)      +M6)/gamma(q6+2);
end
t=0:0.01:20;
figure(1)
plot(t,ww,'Linewidth',1);
xlabel('t/s','fontsize',15);
ylabel('x_4','fontsize',15);
figure(2)
plot(t,gg,'Linewidth',1); 
xlabel('t/s','fontsize',15);
ylabel('x_5','fontsize',15);
figure(3)
plot(t,mm,'Linewidth',1); 
xlabel('t/s','fontsize',15);
ylabel('x_6','fontsize',15);
end
