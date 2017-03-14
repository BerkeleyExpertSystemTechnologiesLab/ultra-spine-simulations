function BodyForceReader(A,q)
bodies = size(A,1)/6;

cableOuts = A*q;
disp(' ')
for a = 1:bodies
    bodyOuts = cableOuts(6*a-5:6*a);
    disp(['For body ' num2str(a) ' cable forces produce'])
    disp(['X Force:' num2str(-bodyOuts(1))])
    disp(['Y Force:' num2str(-bodyOuts(2))])
    disp(['Z Force:' num2str(-bodyOuts(3))])
    disp(['X Moment:' num2str(bodyOuts(4))])
    disp(['Y Moment:' num2str(bodyOuts(5))])
    disp(['Z Moment:' num2str(bodyOuts(6))])
    disp(' ')
end



end