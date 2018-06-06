function d2xi_solved_un = two_d_spine_accel(in1,in2)
%TWO_D_SPINE_ACCEL
%    D2XI_SOLVED_UN = TWO_D_SPINE_ACCEL(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    05-Jun-2018 17:36:31

tensions_un1 = in2(1,:);
tensions_un2 = in2(2,:);
tensions_un3 = in2(3,:);
tensions_un4 = in2(4,:);
xi1 = in1(1,:);
xi2 = in1(2,:);
xi3 = in1(3,:);
xi6 = in1(6,:);
t2 = cos(xi3);
t3 = sin(xi3);
t4 = sqrt(3.0);
t5 = xi6.^2;
t6 = t3.^2;
t7 = pi.*(1.0./6.0);
t8 = t7+xi3;
t9 = cos(t8);
t10 = pi.*(1.0./3.0);
t11 = t10+xi3;
t12 = sin(t11);
t13 = t2.^2;
t14 = -t10+xi3;
t15 = cos(t14);
t16 = cos(t11);
t17 = t6+t13-3.6e1;
t18 = 1.0./t17;
d2xi_solved_un = [t18.*(t3.*t5.*1.836e3-t4.*tensions_un1.*5.4e4+t4.*tensions_un2.*5.4e4+t9.*tensions_un1.*1.08e5+t9.*tensions_un3.*1.08e5-t12.*tensions_un2.*1.08e5-t12.*tensions_un4.*1.08e5-tensions_un1.*xi1.*7.2e5-tensions_un2.*xi1.*7.2e5-tensions_un3.*xi1.*7.2e5-tensions_un4.*xi1.*7.2e5-t3.*t5.*t6.*5.1e1-t3.*t5.*t13.*5.1e1+t2.*t3.*tensions_un1.*2.25e4+t2.*t3.*tensions_un2.*2.25e4+t2.*t3.*tensions_un3.*1.5e3+t2.*t3.*tensions_un4.*1.5e3+t4.*t6.*tensions_un1.*1.5e3-t4.*t6.*tensions_un2.*1.5e3+t2.*t9.*tensions_un3.*1.2e4-t6.*t9.*tensions_un1.*3.0e3-t2.*t12.*tensions_un4.*1.2e4-t6.*t9.*tensions_un3.*3.0e3+t6.*t12.*tensions_un2.*3.0e3+t6.*t12.*tensions_un4.*3.0e3+t6.*tensions_un1.*xi1.*2.0e4+t6.*tensions_un2.*xi1.*2.0e4+t6.*tensions_un3.*xi1.*2.0e4+t6.*tensions_un4.*xi1.*2.0e4+t13.*tensions_un1.*xi1.*8.0e4+t13.*tensions_un2.*xi1.*8.0e4+t2.*t3.*t15.*tensions_un1.*3.0e3+t2.*t3.*t15.*tensions_un3.*3.0e3+t2.*t3.*t16.*tensions_un2.*3.0e3+t2.*t3.*t16.*tensions_un4.*3.0e3+t2.*t3.*tensions_un1.*xi2.*6.0e4+t2.*t3.*tensions_un2.*xi2.*6.0e4-t2.*t3.*tensions_un3.*xi2.*2.0e4-t2.*t3.*tensions_un4.*xi2.*2.0e4-t2.*t9.*tensions_un3.*xi2.*1.6e5+t2.*t12.*tensions_un4.*xi2.*1.6e5-t4.*t13.*tensions_un1.*xi2.*8.0e4+t2.*t15.*tensions_un3.*xi1.*1.6e5+t4.*t13.*tensions_un2.*xi2.*8.0e4+t2.*t16.*tensions_un4.*xi1.*1.6e5+t2.*t3.*t4.*tensions_un1.*xi1.*8.0e4-t2.*t3.*t4.*tensions_un2.*xi1.*8.0e4).*(-3.676470588235294e-4);t18.*(t6.*2.6656e4+t13.*2.6656e4-tensions_un1.*5.4e4-tensions_un2.*5.4e4+tensions_un3.*5.4e4+tensions_un4.*5.4e4-t2.*t5.*1.836e3+t6.*tensions_un1.*2.4e4+t6.*tensions_un2.*2.4e4+t13.*tensions_un1.*1.5e3+t13.*tensions_un2.*1.5e3-t13.*tensions_un3.*1.5e3+t15.*tensions_un1.*1.08e5-t13.*tensions_un4.*1.5e3+t15.*tensions_un3.*1.08e5+t16.*tensions_un2.*1.08e5+t16.*tensions_un4.*1.08e5-tensions_un1.*xi2.*7.2e5-tensions_un2.*xi2.*7.2e5-tensions_un3.*xi2.*7.2e5-tensions_un4.*xi2.*7.2e5+t2.*t5.*t6.*5.1e1+t2.*t5.*t13.*5.1e1+t3.*t9.*tensions_un3.*1.2e4-t3.*t12.*tensions_un4.*1.2e4-t13.*t15.*tensions_un1.*3.0e3-t13.*t15.*tensions_un3.*3.0e3-t13.*t16.*tensions_un2.*3.0e3-t13.*t16.*tensions_un4.*3.0e3+t6.*tensions_un1.*xi2.*8.0e4+t6.*tensions_un2.*xi2.*8.0e4+t13.*tensions_un1.*xi2.*2.0e4+t13.*tensions_un2.*xi2.*2.0e4+t13.*tensions_un3.*xi2.*2.0e4+t13.*tensions_un4.*xi2.*2.0e4-t2.*t3.*t4.*tensions_un1.*1.5e3+t2.*t3.*t4.*tensions_un2.*1.5e3+t2.*t3.*t9.*tensions_un1.*3.0e3+t2.*t3.*t9.*tensions_un3.*3.0e3-t2.*t3.*t12.*tensions_un2.*3.0e3-t2.*t3.*t12.*tensions_un4.*3.0e3+t2.*t3.*tensions_un1.*xi1.*6.0e4+t2.*t3.*tensions_un2.*xi1.*6.0e4-t2.*t3.*tensions_un3.*xi1.*2.0e4-t2.*t3.*tensions_un4.*xi1.*2.0e4+t4.*t6.*tensions_un1.*xi1.*8.0e4-t4.*t6.*tensions_un2.*xi1.*8.0e4-t3.*t9.*tensions_un3.*xi2.*1.6e5+t3.*t12.*tensions_un4.*xi2.*1.6e5+t3.*t15.*tensions_un3.*xi1.*1.6e5+t3.*t16.*tensions_un4.*xi1.*1.6e5-t2.*t3.*t4.*tensions_un1.*xi2.*8.0e4+t2.*t3.*t4.*tensions_un2.*xi2.*8.0e4-9.59616e5).*(-3.676470588235294e-4);t18.*(t3.*tensions_un1.*4.5e1+t3.*tensions_un2.*4.5e1+t3.*tensions_un3.*3.0+t3.*tensions_un4.*3.0+t9.*tensions_un3.*2.4e1-t12.*tensions_un4.*2.4e1-t2.*t4.*tensions_un1.*3.0+t2.*t4.*tensions_un2.*3.0+t2.*t9.*tensions_un1.*6.0+t2.*t9.*tensions_un3.*6.0-t2.*t12.*tensions_un2.*6.0-t2.*t12.*tensions_un4.*6.0+t3.*t15.*tensions_un1.*6.0+t3.*t15.*tensions_un3.*6.0+t3.*t16.*tensions_un2.*6.0+t3.*t16.*tensions_un4.*6.0+t2.*tensions_un1.*xi1.*1.2e2+t2.*tensions_un2.*xi1.*1.2e2-t2.*tensions_un3.*xi1.*4.0e1+t3.*tensions_un1.*xi2.*1.2e2-t2.*tensions_un4.*xi1.*4.0e1+t3.*tensions_un2.*xi2.*1.2e2-t3.*tensions_un3.*xi2.*4.0e1-t3.*tensions_un4.*xi2.*4.0e1-t9.*tensions_un3.*xi2.*3.2e2+t12.*tensions_un4.*xi2.*3.2e2+t15.*tensions_un3.*xi1.*3.2e2+t16.*tensions_un4.*xi1.*3.2e2-t2.*t4.*tensions_un1.*xi2.*1.6e2+t3.*t4.*tensions_un1.*xi1.*1.6e2+t2.*t4.*tensions_un2.*xi2.*1.6e2-t3.*t4.*tensions_un2.*xi1.*1.6e2).*(5.0e2./5.1e1)];
