function lengths = two_d_spine_lengths(in1)
%TWO_D_SPINE_LENGTHS
%    LENGTHS = TWO_D_SPINE_LENGTHS(IN1)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    10-Dec-2016 14:59:59

xi1 = in1(1,:);
xi2 = in1(2,:);
xi3 = in1(3,:);
t4 = xi1.*4.0e1;
t6 = sqrt(3.0);
t7 = t6.*3.0;
t12 = pi.*(1.0./6.0);
t13 = t12+xi3;
t14 = cos(t13);
t2 = t4+t7-t14.*6.0;
t5 = pi.*(1.0./3.0);
t9 = xi2.*4.0e1;
t16 = -t5+xi3;
t17 = cos(t16);
t18 = t17.*6.0;
t3 = t9-t18+3.0;
t10 = t5+xi3;
t24 = sin(t10);
t8 = t4-t7+t24.*6.0;
t20 = cos(t10);
t21 = t20.*6.0;
t11 = t9-t21+3.0;
t23 = xi1.*2.0e1;
t15 = t14.*3.0-t23;
t19 = -t9+t18+3.0;
t22 = -t9+t21+3.0;
t25 = t23+t24.*3.0;
lengths = [sqrt(t2.^2+t3.^2).*(1.0./4.0e1);sqrt(t8.^2+t11.^2).*(1.0./4.0e1);sqrt(t15.^2.*4.0+t19.^2).*(1.0./4.0e1);sqrt(t22.^2+t25.^2.*4.0).*(1.0./4.0e1)];
