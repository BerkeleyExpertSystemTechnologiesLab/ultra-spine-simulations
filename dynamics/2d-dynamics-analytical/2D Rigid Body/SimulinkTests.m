%%Simulink Test
clear
clc
close all

tstep = 5;
m = 3;
% m = 10000;
% I = (0.5^2*0.33)*3;
I = 10000;
l = 1;
r0 = [-1 0.5];
k =  1500;
ks = 1500;

A = Spine2(1,109.5);
x2 = [0 2];
x1 = [-1 0 0];
l1 = [-1 A.bars(1).lx(2), A.bars(1).ly(2)];
l2 = [-1 A.bars(2).lx(2), A.bars(2).ly(2)];
l3 = [-1 A.bars(3).lx(2), A.bars(3).ly(2)];

d = 0.0257;

theta0 = 0;
theta = [-1 A.th];

simout = sim('Spine2sim');


%%
figure;
subplot(3,1,1)
plot(pEnergy.Time,pEnergy.Data + sum(kEnergy.Data,2))
title('Total Energy vs Time: No Moment')

subplot(3,1,2)
plot(pEnergy.Time,pEnergy.Data)
title('Potential Energy vs Time')

subplot(3,1,3)
plot(pEnergy.Time,sum(kEnergy.Data,2))
title('Total Kinetic Energy vs Time');

%%
figure;

numz=0;

for i = 2:length(Xmain.Data)
%     pause(Xmain.Time(i)-Xmain.Time(i-1));
    L(1).x = [x1(2),a1.Data(end,1)];
    L(2).x = [x1(2),a2.Data(end,1)];
    L(3).x = [x1(2),a3.Data(end,1)];
    L(1).y = [x1(3),a1.Data(end,2)];
    L(2).y = [x1(3),a2.Data(end,2)];
    L(3).y = [x1(3),a3.Data(end,2)];
    L(4).x = [Xmain.Data(i,1),b1.Data(i,1)];
    L(5).x = [Xmain.Data(i,1),b2.Data(i,1)];
    L(6).x = [Xmain.Data(i,1),b3.Data(i,1)];
    L(4).y = [Xmain.Data(i,2),b1.Data(i,2)];
    L(5).y = [Xmain.Data(i,2),b2.Data(i,2)];
    L(6).y = [Xmain.Data(i,2),b3.Data(i,2)];
    
    if round(Xmain.Time(i),3) ~= round(Xmain.Time(i-1),3)
        for n = 1:3
            plot(L(n).x,L(n).y,'b')
            hold on
        end
        
        for n = 4:6
            plot(L(n).x,L(n).y,'r');
        end
        axis([-4 4 -4 4]);
        hold off
        legend(num2str(Xmain.Time(i)))
        numz = numz+1;
        drawnow;
        Frames(numz) = getframe;
    end
end

%%
i = 387;
t = 2;
while t < 3
    i = i + 1;
    t = Xmain.Time(i);
end

%%

% for n = 1:length(xSaved.data)
%     
% end

% %%
% xSaved = [];
% LL1 = []; LL2 = []; LL3 = [];
% for n=1:1000
%     simout = sim('Spine2sim');
%     x2 = Xmain.Data(end,:);
%     theta0 = THTA.Data(end);
%     xSaved = [xSaved;x2];
%     
%     LL1 = [LL1;b1.Data(end,:)];
%     LL2 = [LL2;b2.Data(end,:)];
%     LL3 = [LL3;b3.Data(end,:)];
% end
% 
% %%
% 
% for i = 1:1000
%     pause(0.05);
%     L(1).x = [x1(2),a1.Data(end,1)];
%     L(2).x = [x1(2),a2.Data(end,1)];
%     L(3).x = [x1(2),a3.Data(end,1)];
%     L(1).y = [x1(3),a1.Data(end,2)];
%     L(2).y = [x1(3),a2.Data(end,2)];
%     L(3).y = [x1(3),a3.Data(end,2)];
%     L(4).x = [xSaved(i,1),LL1(i,1)];
%     L(5).x = [xSaved(i,1),LL2(i,1)];
%     L(6).x = [xSaved(i,1),LL3(i,1)];
%     L(4).y = [xSaved(i,2),LL1(i,2)];
%     L(5).y = [xSaved(i,2),LL2(i,2)];
%     L(6).y = [xSaved(i,2),LL3(i,2)];
%     
%     
%     
%     for n = 1:3
%         plot(L(n).x,L(n).y,'b')
%         hold on
%     end
%     
%     for n = 4:6
%         plot(L(n).x,L(n).y,'r');
%     end
%     
%     hold off
%     
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% L(1).x = [x1(2),a1.Data(end,1)];
% L(2).x = [x1(2),a2.Data(end,1)];
% L(3).x = [x1(2),a3.Data(end,1)];
% L(1).y = [x1(3),a1.Data(end,2)];
% L(2).y = [x1(3),a2.Data(end,2)];
% L(3).y = [x1(3),a3.Data(end,2)];
% L(4).x = [Xmain.Data(end,1),b1.Data(end,1)];
% L(5).x = [Xmain.Data(end,1),b2.Data(end,1)];
% L(6).x = [Xmain.Data(end,1),b3.Data(end,1)];
% L(4).y = [Xmain.Data(end,2),b1.Data(end,2)];
% L(5).y = [Xmain.Data(end,2),b2.Data(end,2)];
% L(6).y = [Xmain.Data(end,2),b3.Data(end,2)];
% 
% 
% hold on
% for n = 1:3
%     plot(L(n).x,L(n).y,'b')
% end
% 
% for n = 4:6
%     plot(L(n).x,L(n).y,'r');
% end
