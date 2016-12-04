function dlengths_dt = two_d_spine_dlengths_dt(x2,z2,T2,dx2,dz2,dT2)
%TWO_D_SPINE_DLENGTHS_DT
%    DLENGTHS_DT = TWO_D_SPINE_DLENGTHS_DT(X2,Z2,T2,DX2,DZ2,DT2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    27-Nov-2016 16:36:17

t4 = pi.*(1.0./6.0);
t5 = T2+t4;
t6 = cos(t5);
t7 = t6.*(3.0./2.0e1);
t8 = sqrt(3.0);
t9 = t8.*(3.0./4.0e1);
t10 = -t7+t9+x2;
t2 = abs(t10);
t11 = pi.*(1.0./3.0);
t12 = T2-t11;
t13 = cos(t12);
t14 = t13.*(3.0./2.0e1);
t15 = -t14+z2+3.0./4.0e1;
t3 = abs(t15);
t16 = sign(t15);
t17 = sign(t10);
t18 = sin(T2);
t19 = cos(T2);
t22 = T2-t4;
t23 = cos(t22);
t24 = t23.*(3.0./2.0e1);
t25 = -t9+t24+x2;
t20 = abs(t25);
t26 = T2+t11;
t27 = cos(t26);
t28 = t27.*(3.0./2.0e1);
t29 = -t28+z2+3.0./4.0e1;
t21 = abs(t29);
t30 = sign(t29);
t31 = sign(t25);
t34 = -t7+x2;
t32 = abs(t34);
t35 = -t14+z2-3.0./4.0e1;
t33 = abs(t35);
t36 = sign(t34);
t37 = sign(t35);
t40 = t24+x2;
t38 = abs(t40);
t41 = -t28+z2-3.0./4.0e1;
t39 = abs(t41);
t42 = sign(t40);
t43 = sign(t41);
dlengths_dt = reshape([0.0,0.0,0.0,0.0,1.0./sqrt(t2.^2+t3.^2).*(dx2.*t2.*t17.*4.0e1+dz2.*t3.*t16.*4.0e1+dT2.*t3.*t16.*t18.*3.0+dT2.*t2.*t17.*t19.*3.0+dT2.*t2.*t8.*t17.*t18.*3.0-dT2.*t3.*t8.*t16.*t19.*3.0).*(1.0./4.0e1),1.0./sqrt(t20.^2+t21.^2).*(dx2.*t20.*t31.*4.0e1+dz2.*t21.*t30.*4.0e1+dT2.*t18.*t21.*t30.*3.0+dT2.*t19.*t20.*t31.*3.0-dT2.*t8.*t18.*t20.*t31.*3.0+dT2.*t8.*t19.*t21.*t30.*3.0).*(1.0./4.0e1),1.0./sqrt(t32.^2+t33.^2).*(dx2.*t32.*t36.*4.0e1+dz2.*t33.*t37.*4.0e1+dT2.*t19.*t32.*t36.*3.0+dT2.*t18.*t33.*t37.*3.0+dT2.*t8.*t18.*t32.*t36.*3.0-dT2.*t8.*t19.*t33.*t37.*3.0).*(1.0./4.0e1),1.0./sqrt(t38.^2+t39.^2).*(dx2.*t38.*t42.*4.0e1+dz2.*t39.*t43.*4.0e1+dT2.*t19.*t38.*t42.*3.0+dT2.*t18.*t39.*t43.*3.0-dT2.*t8.*t18.*t38.*t42.*3.0+dT2.*t8.*t19.*t39.*t43.*3.0).*(1.0./4.0e1)],[4,2]);
