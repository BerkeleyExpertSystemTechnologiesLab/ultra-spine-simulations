clear variables
close all
clc
lim =0.175; %plot limits
del_t = 0.01; %time step for dynamics
nn = 3; %frame divisor
stringEnable =1; %enable string plotting
l = 0.159*2; %meters
h = 0.23; % meters

x=0.0; y=-0.0; z=0.15; T=0.0; G=0.0; P=0.0;
dx=0.01; dy=0; dz=0; dT=0; dG=0; dP=0;

load bulletData
u1 = RestLength0/100;
u2 = RestLength1/100;
u3 = RestLength2/100;
u4 = RestLength3/100;
u5 = RestLength4/100;
u6 = RestLength5/100;
u7 = RestLength6/100;
u8 = RestLength7/100;

Figs = figure('units','normalized','outerposition',[0 0 1 1]);

M = struct('cdata', cell(1,round(length(u1)/10)), 'colormap', cell(1,round(length(u1)/10)));
rad = 0.01;

Tetra2 =  [...
    0   -l/2+0.0048      -h/2 ; %A
    0    l/2-0.0048      -h/2 ; %B
    -l/2+0.0048  0         h/2 ; %C
    l/2-0.0048  0         h/2 ];%D
Tetra1 =  [...
    0   -l/2+0.0048      -h/2 ; %A
    0    l/2-0.0048      -h/2 ; %B
    -l/2+0.0048  0         h/2 ; %C
    l/2-0.0048  0         h/2 ];%D

Tetra =  [...
    0   -l/2-0.005      -h/2 ; %A
    0    l/2+0.005      -h/2 ; %B
    -l/2-0.005  0         h/2 ; %C
    l/2+0.005  0         h/2 ];%D
cmaps = winter(512);
ax = axes();
[T_base,h_base] = plotTetra(Tetra,rad,ax);
[T_top,h_top] = plotTetra(Tetra,rad,ax);
if(stringEnable)
    
    anchor1=[0 0 rad];
    anchor2=[0 0 rad];
    
    String_pts = [...
        (Tetra1(1,:)+anchor1);
        (Tetra2(1,:)-anchor2);
        (Tetra2(1,:));
        (Tetra2(2,:));
        (Tetra2(2,:)-anchor2);
        (Tetra1(2,:)+anchor1);
        (Tetra1(2,:));
        (Tetra1(3,:));
        (Tetra1(3,:)+anchor1);
        (Tetra2(3,:)-anchor2);
        (Tetra2(3,:));
        (Tetra2(4,:));
        (Tetra2(4,:)-anchor2);
        (Tetra1(4,:)+anchor1);
        (Tetra1(4,:)-anchor1);
        (Tetra2(1,:)+anchor2);
        (Tetra1(3,:)-anchor1);
        (Tetra2(2,:)+anchor2);
        (Tetra1(4,:)-anchor1)];
    string_handle=plot3(String_pts(:,1),String_pts(:,2),String_pts(:,3),'LineWidth',2,'Color','m');
    linkdata on
end
grid on;
axis equal;
xlim([-lim-0.085,lim+0.085])
ylim([-lim,lim])
zlim([-0.15,2*lim])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
colormap(cmaps(1:256,:))
shading interp
light
lighting phong
hold on
view(3)
title('Euler-Lagrange Model')
m = 256;  % Number of colors in the current
stateDy=zeros(round(length(u1)/10),12);
time = 0:del_t:del_t*length(u1);
for i=2:(length(u1)/10)
    
    if mod(i,nn) == 0
        tic
        while toc<del_t*nn
        end
        disp(toc)
        
        
        
        RR1 =  getHG_Tform(x,y,z,T,G,P);
        set(T_top,'Matrix',RR1);
        
        if stringEnable
            anchor2 = [0 0 rad 1];
            Tetra2 =  [...
                0   -l/2+0.0048      -h/2  1; %A
                0    l/2-0.0048      -h/2  1; %B
                -l/2+0.0048  0         h/2  1; %C
                l/2-0.0048  0         h/2  1];%D
            Tetra2 = RR1*Tetra2';
            Tetra2 = Tetra2';
            RR1(1:3,4) = [0; 0; 0];
            anchor2 = RR1*anchor2';
            anchor2 = anchor2(1:3)';
            Tetra2 = Tetra2(:,1:3);
            String_pts = [...
                (Tetra1(1,:)+anchor1);
                (Tetra2(1,:)-anchor2);
                (Tetra2(1,:));
                (Tetra2(2,:));
                (Tetra2(2,:)-anchor2);
                (Tetra1(2,:)+anchor1);
                (Tetra1(2,:));
                (Tetra1(3,:));
                (Tetra1(3,:)+anchor1);
                (Tetra2(3,:)-anchor2);
                (Tetra2(3,:));
                (Tetra2(4,:));
                (Tetra2(4,:)-anchor2);
                (Tetra1(4,:)+anchor1);
                (Tetra1(4,:)-anchor1);
                (Tetra2(1,:)+anchor2);
                (Tetra1(3,:)-anchor1);
                (Tetra2(2,:)+anchor2);
                (Tetra1(4,:)-anchor1)];
            refreshdata(string_handle);
        end
        drawnow;
    end
    
    dumX = [x,y,z,T,G,P,dx,dy,dz,dT,dG,dP];
    Te = getTensions(dumX(1),dumX(2),dumX(3),dumX(4),dumX(5),dumX(6),dumX(7),dumX(8),dumX(9),dumX(10),dumX(11),dumX(12),u1(i*10),u2(i*10),u3(i*10),u4(i*10),u5(i*10),u6(i*10),u7(i*10),u8(i*10));
    K1(1:6) = dumX(7:12);
    K1(7:12) = duct_accel(dumX(1),dumX(2),dumX(3),dumX(4),dumX(5),dumX(6),dumX(7),dumX(8),dumX(9),dumX(10),dumX(11),dumX(12),Te(1),Te(2),Te(3),Te(4),Te(5),Te(6),Te(7),Te(8));
    Te = getTensions(dumX(1),dumX(2),dumX(3),dumX(4),dumX(5),dumX(6),dumX(7),dumX(8),dumX(9),dumX(10),dumX(11),dumX(12),u1(i*10),u2(i*10),u3(i*10),u4(i*10),u5(i*10),u6(i*10),u7(i*10),u8(i*10));
    dumX = [x,y,z,T,G,P,dx,dy,dz,dT,dG,dP]+K1*del_t/2;
    K2(1:6) =  dumX(7:12);
    K2(7:12) =  duct_accel(dumX(1),dumX(2),dumX(3),dumX(4),dumX(5),dumX(6),dumX(7),dumX(8),dumX(9),dumX(10),dumX(11),dumX(12),Te(1),Te(2),Te(3),Te(4),Te(5),Te(6),Te(7),Te(8));
    Te = getTensions(dumX(1),dumX(2),dumX(3),dumX(4),dumX(5),dumX(6),dumX(7),dumX(8),dumX(9),dumX(10),dumX(11),dumX(12),u1(i*10),u2(i*10),u3(i*10),u4(i*10),u5(i*10),u6(i*10),u7(i*10),u8(i*10));
    dumX = [x,y,z,T,G,P,dx,dy,dz,dT,dG,dP]+K2*del_t/2;
    K3(1:6) =  dumX(7:12);
    K3(7:12) =  duct_accel(dumX(1),dumX(2),dumX(3),dumX(4),dumX(5),dumX(6),dumX(7),dumX(8),dumX(9),dumX(10),dumX(11),dumX(12),Te(1),Te(2),Te(3),Te(4),Te(5),Te(6),Te(7),Te(8));
    Te = getTensions(dumX(1),dumX(2),dumX(3),dumX(4),dumX(5),dumX(6),dumX(7),dumX(8),dumX(9),dumX(10),dumX(11),dumX(12),u1(i*10),u2(i*10),u3(i*10),u4(i*10),u5(i*10),u6(i*10),u7(i*10),u8(i*10));
    dumX = [x,y,z,T,G,P,dx,dy,dz,dT,dG,dP]+K3*del_t;
    K4(1:6) =  dumX(7:12);
    K4(7:12) =  duct_accel(dumX(1),dumX(2),dumX(3),dumX(4),dumX(5),dumX(6),dumX(7),dumX(8),dumX(9),dumX(10),dumX(11),dumX(12),Te(1),Te(2),Te(3),Te(4),Te(5),Te(6),Te(7),Te(8));
    stateDy(i,:) = [x,y,z,T,G,P,dx,dy,dz,dT,dG,dP]+del_t/6*(K1+2*K2+2*K3+K4);
    
    x = stateDy(i,1);
    y = stateDy(i,2);
    z = stateDy(i,3);
    T = stateDy(i,4);
    G = stateDy(i,5);
    P = stateDy(i,6);
    dx = stateDy(i,7);
    dy = stateDy(i,8);
    dz = stateDy(i,9);
    dT = stateDy(i,10);
    dG = stateDy(i,11);
    dP = stateDy(i,12);
end



