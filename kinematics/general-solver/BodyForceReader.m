function BodyForceReader(q,A,p)
% Takes equilibrium matrix and cable force densities, prints out forces and
% moments from cables on bodies.
% Specifically, this compares the resulting forces/moments on the bodies
% due to the cable tensions q against the desired applied forces, p.

bodies = size(A,1)/6;

cableOuts = A*q;
disp(' ')
disp('Comparing the forces/moments due to cable tensions against the desired forces...');
for a = 1:bodies
    indices = a:bodies:a+5*bodies;
    bodyOuts = cableOuts(indices);
    disp(['For body ' num2str(a) ' cable forces produce'])
    disp(['X Force:' num2str(-bodyOuts(1)) ', compare to the desired equilibrium force:' num2str(-p(indices(1)))])
    disp(['Y Force:' num2str(-bodyOuts(2)) ', compare to the desired equilibrium force:' num2str(-p(indices(2)))])
    disp(['Z Force:' num2str(-bodyOuts(3)) ', compare to the desired equilibrium force:' num2str(-p(indices(3)))])
    disp(['X Moment:' num2str(bodyOuts(4)) ', compare to the desired equilibrium force:' num2str(p(indices(4)))])
    disp(['Y Moment:' num2str(bodyOuts(5)) ', compare to the desired equilibrium force:' num2str(p(indices(5)))])
    disp(['Z Moment:' num2str(bodyOuts(6)) ', compare to the desired equilibrium force:' num2str(p(indices(6)))])
    disp(' ')
end

disp('If any of these are different, the bodies will not stay still.');
disp(' ');

end