function tensions_sub = two_d_spine_tensions(in1,in2)
%TWO_D_SPINE_TENSIONS
%    TENSIONS_SUB = TWO_D_SPINE_TENSIONS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    04-Dec-2016 16:37:32

u1 = in2(1,:);
u2 = in2(2,:);
u3 = in2(3,:);
u4 = in2(4,:);
xi1 = in1(1,:);
xi2 = in1(2,:);
xi3 = in1(3,:);
xi4 = in1(4,:);
xi5 = in1(5,:);
xi6 = in1(6,:);
t4 = xi1.*4.0e1;
t5 = pi.*(1.0./6.0);
t6 = t5+xi3;
t7 = cos(t6);
t8 = t7.*6.0;
t9 = sqrt(3.0);
t10 = t9.*3.0;
t2 = t4-t8+t10;
t12 = xi2.*4.0e1;
t13 = pi.*(1.0./3.0);
t14 = -t13+xi3;
t15 = cos(t14);
t16 = t15.*6.0;
t3 = t12-t16+3.0;
t11 = t2.^2;
t17 = t3.^2;
t18 = t11+t17;
t19 = sin(xi3);
t20 = cos(xi3);
t22 = t13+xi3;
t24 = sin(t22);
t25 = t24.*6.0;
t21 = t4-t10+t25;
t27 = cos(t22);
t28 = t27.*6.0;
t23 = t12-t28+3.0;
t26 = t21.^2;
t29 = t23.^2;
t30 = t26+t29;
t31 = xi5.*3.0e1;
t32 = t9.*xi4.*3.0e1;
t33 = xi1.*xi4.*4.0e2;
t34 = xi2.*xi5.*4.0e2;
t35 = t19.*xi4.*3.0e1;
t36 = t19.*xi6.*9.0;
t37 = t20.*xi1.*xi6.*3.0e1;
t38 = t19.*xi2.*xi6.*3.0e1;
t39 = t9.*t19.*xi1.*xi6.*3.0e1;
t42 = xi1.*2.0e1;
t43 = t7.*3.0;
t40 = -t42+t43;
t41 = -t12+t16+3.0;
t44 = t42-t43;
t45 = t41.^2;
t46 = -t12+t28+3.0;
t51 = t24.*3.0;
t47 = t42+t51;
t48 = xi1.*xi4.*8.0e2;
t49 = xi2.*xi5.*8.0e2;
t50 = t46.^2;
t52 = t47.^2;
t53 = t52.*4.0;
t54 = t50+t53;
tensions_sub = [u1.*-2.0e3+1.0./sqrt(t18).*(t31+t32+t33+t34+t35+t36+t37+t38+t39-t20.*xi5.*3.0e1-t9.*t19.*xi5.*3.0e1-t9.*t20.*xi4.*3.0e1-t9.*t20.*xi2.*xi6.*3.0e1).*(5.0./2.0)+sqrt(t18).*5.0e1;u2.*-2.0e3+1.0./sqrt(t30).*(t31-t32+t33+t34+t35+t36+t37+t38-t39-t20.*xi5.*3.0e1+t9.*t19.*xi5.*3.0e1+t9.*t20.*xi4.*3.0e1+t9.*t20.*xi2.*xi6.*3.0e1).*(5.0./2.0)+sqrt(t30).*5.0e1;u3.*-2.0e3+1.0./sqrt(t45+t40.^2.*4.0).*(t48+t49-xi5.*6.0e1+t19.*xi4.*6.0e1-t19.*xi6.*(9.0./2.0)-t20.*xi5.*6.0e1-t9.*t19.*xi5.*6.0e1-t9.*t20.*xi4.*6.0e1+t9.*t20.*xi6.*(9.0./2.0)+t19.*xi2.*xi6.*6.0e1+t20.*xi1.*xi6.*6.0e1+t9.*t19.*xi1.*xi6.*6.0e1-t9.*t20.*xi2.*xi6.*6.0e1).*(5.0./4.0)+sqrt(t45+t44.^2.*4.0).*5.0e1;u4.*-2.0e3+1.0./sqrt(t54).*(t48+t49-xi5.*6.0e1+t24.*xi4.*1.2e2-t24.*xi6.*9.0-t27.*xi5.*1.2e2+t24.*xi2.*xi6.*1.2e2+t27.*xi1.*xi6.*1.2e2).*(5.0./4.0)+sqrt(t54).*5.0e1];
