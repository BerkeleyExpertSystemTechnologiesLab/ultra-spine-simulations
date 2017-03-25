generalInverseKinematicsSolverv4
load('verifiedOutput.mat')
error = norm(A-ACheck)+norm(p-pCheck)+norm(tensions-tensionsCheck);

Fx = -q(:).*dx(allCables);
Fy = -q(:).*dy(allCables);
Fz = -q(:).*dz(allCables);
sum(Fx(1:8))
sum(Fy(1:8))
sum(Fz(1:8))
