clear
load('datav2.mat')

plotter = @(a,b,c) line([x(a) x(b)],[y(a) y(b)],[z(a) z(b)],'color',c);

for a = 1:32
    indices = find(C(a,:));
   plotter(indices(1),indices(2),[0,0,1])
end
   

for a = 33:39
    indices = find(C(a,:));
   plotter(indices(1),indices(2),[1,0,0])
end

comx = zeros(7,1);
comy = zeros(7,1);
comz = zeros(7,1);

for a = 1:7
    comx(a) = (x(2*a-1) + x(2*a)) /2;
    comy(a) = (y(2*a-1) + y(2*a)) /2;
    comz(a) = (z(2*a-1) + z(2*a)) /2;
end

line(comx,comy,comz,'Marker','o','linestyle','none')

x = [x;comx];
y = [y;comy];
z = [z;comz];

axis equal