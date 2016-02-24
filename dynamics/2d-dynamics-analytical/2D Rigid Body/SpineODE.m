function dout = SpineODE(t,xx,l1,l2,l3,xl,a,k,I,m,r0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = xx(1);
y = xx(3);
theta = xx(5);

L{1}(1) = x + xl.*cos(a(1) - theta);
L{2}(1) = x + xl.*cos(a(2) - theta);
L{3}(1) = x + xl.*cos(a(3) - theta);

L{1}(2) = y + xl.*sin(a(1) - theta);
L{2}(2) = y + xl.*sin(a(2) - theta);
L{3}(2) = y + xl.*sin(a(3) - theta);

r1 = L{1} - l1;
r2 = L{2} - l2;
r3 = L{3} - l3;
r4 = L{3} - l1;
r5 = L{3} - l2;

u1 = r1./norm(r1);
u2 = r2./norm(r2);
u3 = r3./norm(r3);
u4 = r4./norm(r4);
u5 = r5./norm(r5);

%% Force Calculation
for i = 1:2
    F{1}(i) = k*(norm(r1) - r0).*u1(i);
    F{2}(i) = k*(norm(r2) - r0).*u2(i);
    F{3}(i) = k*(norm(r3) - r0).*u3(i);
    F{4}(i) = k*(norm(r4) - r0).*u4(i);
    F{5}(i) = k*(norm(r5) - r0).*u5(i);
end

%Sum of the forces
Fsum = 0;
for n = 1:5
    Fsum = Fsum + F{n};
end

ax = -Fsum(1)./m;
ax = round(ax,8);
ay = -Fsum(2)./m;
ay = round(ay,12);


%% Moment Calculation
for n = 1:5
   F1{n} = [F{n} 0]; 
end

for n = 1:3
    L1{n} = [L{n} 0];
    
%     M1 = cross(L1{n},F1{n});
    M1 = (L1{n}(1)-x(1))*F1{n}(2) - (L1{n}(2)-y(1))*F1{n}(1);
%     keyboard
    M(n) = M1(end);
end

M1 = cross(L1{3},F1{4});
M(4) = M1(end);
M1 = cross(L1{3},F1{5});
M(5) = M1(end);

Msum = sum(M);
Msum = round(Msum,8);
alpha = Msum./I;

%% Build dout

dout(1) = xx(2);
dout(2) = ax;
dout(3) = xx(4);
dout(4) = ay;
dout(5) = xx(6);
dout(6) = alpha;

dout = dout';

% keyboard
