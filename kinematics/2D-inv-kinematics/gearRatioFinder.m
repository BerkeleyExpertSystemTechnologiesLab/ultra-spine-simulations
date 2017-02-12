% ChanWoo

% Finding gear ratio

% stringLengthHistoryVert(A,B,C)
% A: Connection number
% B: Cable Number
% C: Time Frame

%% Analyze changes in cable length over entire time frame

%VertCable'i'1 = Longest Cable
%VertCable'i'1 = Second Longest Cable
%VertCable'i'1 = Third Longest Cable
%VertCable'i'1 = Fourth Longest Cable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          
% 0-------------0-------------0-------------0-------------0     vertCablei1 
% 0-------------0-------------0-------------0                   vertCablei2 
% 0-------------0-------------0                                 vertCablei3 
% 0-------------0                                               vertCablei4 
% 
% 0 represents the joint of tetrahedron. (Connection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertCable11 = stringLengthHistoryVert(1,1,:)+stringLengthHistoryVert(2,1,:)+...
              stringLengthHistoryVert(3,1,:)+stringLengthHistoryVert(4,1,:);
vertCable12 = stringLengthHistoryVert(1,1,:)+stringLengthHistoryVert(2,1,:)+...
              stringLengthHistoryVert(3,1,:);
vertCable13 = stringLengthHistoryVert(1,1,:)+stringLengthHistoryVert(2,1,:);
vertCable14 = stringLengthHistoryVert(1,1,:);

vertCable21 = stringLengthHistoryVert(1,2,:)+stringLengthHistoryVert(2,2,:)+...
              stringLengthHistoryVert(3,2,:)+stringLengthHistoryVert(4,2,:);
vertCable22 = stringLengthHistoryVert(1,2,:)+stringLengthHistoryVert(2,2,:)+...
              stringLengthHistoryVert(3,2,:);
vertCable23 = stringLengthHistoryVert(1,2,:)+stringLengthHistoryVert(2,2,:);
vertCable24 = stringLengthHistoryVert(1,2,:);

vertCable31 = stringLengthHistoryVert(1,3,:)+stringLengthHistoryVert(2,3,:)+...
              stringLengthHistoryVert(3,3,:)+stringLengthHistoryVert(4,3,:);
vertCable32 = stringLengthHistoryVert(1,3,:)+stringLengthHistoryVert(2,3,:)+...
              stringLengthHistoryVert(3,3,:);
vertCable33 = stringLengthHistoryVert(1,3,:)+stringLengthHistoryVert(2,3,:);
vertCable34 = stringLengthHistoryVert(1,3,:);

vertCable41 = stringLengthHistoryVert(1,4,:)+stringLengthHistoryVert(2,4,:)+...
              stringLengthHistoryVert(3,4,:)+stringLengthHistoryVert(4,4,:);
vertCable42 = stringLengthHistoryVert(1,4,:)+stringLengthHistoryVert(2,4,:)+...
              stringLengthHistoryVert(3,4,:);
vertCable43 = stringLengthHistoryVert(1,4,:)+stringLengthHistoryVert(2,4,:);
vertCable44 = stringLengthHistoryVert(1,4,:);

%% Analyze in change of length at each time point (delta length / delta time)
for i = 1:num_frames_to_render-1
    dVertCable31(i) = vertCable31(:,:,i+1)-vertCable31(:,:,i);
    dVertCable32(i) = vertCable32(:,:,i+1)-vertCable32(:,:,i);
    dVertCable33(i) = vertCable33(:,:,i+1)-vertCable33(:,:,i);
    dVertCable34(i) = vertCable34(:,:,i+1)-vertCable34(:,:,i);
    
    dVertCable41(i) = vertCable41(:,:,i+1)-vertCable41(:,:,i);
    dVertCable42(i) = vertCable42(:,:,i+1)-vertCable42(:,:,i);
    dVertCable43(i) = vertCable43(:,:,i+1)-vertCable43(:,:,i);
    dVertCable44(i) = vertCable44(:,:,i+1)-vertCable44(:,:,i);
end
%% Plot
%% Constrained
figure(20)
suptitle('Constrained Model')
subplot(2,2,1)
hold on
plot(squeeze(vertCable11(:,:,:)))
plot(squeeze(vertCable12(:,:,:)))
plot(squeeze(vertCable13(:,:,:)))
plot(squeeze(vertCable14(:,:,:)))
title('Changes in Vertical Cable1 Length')
grid on
subplot(2,2,2)
hold on
plot(squeeze(vertCable21(:,:,:)))
plot(squeeze(vertCable22(:,:,:)))
plot(squeeze(vertCable23(:,:,:)))
plot(squeeze(vertCable24(:,:,:)))
title('Changes in Vertical Cable2 Length 2')
grid on
subplot(2,2,3)
hold on
plot(squeeze(vertCable31(:,:,:)))
plot(squeeze(vertCable32(:,:,:)))
plot(squeeze(vertCable33(:,:,:)))
plot(squeeze(vertCable34(:,:,:)))
title('Changes in Vertical Cable3 Length')
grid on
subplot(2,2,4)
hold on
plot(squeeze(vertCable41(:,:,:)))
plot(squeeze(vertCable42(:,:,:)))
plot(squeeze(vertCable43(:,:,:)))
plot(squeeze(vertCable44(:,:,:)))
title('Changes in Vertical Cable4 Length 2')
grid on

% subplot(2,2,3)
% hold on
% plot(squeeze(dVertCable31(:,:,:)))
% plot(squeeze(dVertCable32(:,:,:)))
% plot(squeeze(dVertCable33(:,:,:)))
% plot(squeeze(dVertCable34(:,:,:)))
% title('delta Length of Vertical Cable3')
% grid on
% subplot(2,2,4)
% hold on
% plot(squeeze(dVertCable41(:,:,:)))
% plot(squeeze(dVertCable42(:,:,:)))
% plot(squeeze(dVertCable43(:,:,:)))
% plot(squeeze(dVertCable44(:,:,:)))
% title('delta Length of Vertical Cable4')
% grid on

%% Unconstrained
figure(21)
suptitle('Unconstrained Model')
subplot(2,2,1)
hold on
plot(squeeze(stringLengthHistoryVert(1,1,:)))
plot(squeeze(stringLengthHistoryVert(2,1,:)))
plot(squeeze(stringLengthHistoryVert(3,1,:)))
plot(squeeze(stringLengthHistoryVert(4,1,:)))
grid on
title('Vertical Cable 1')
subplot(2,2,2)
hold on
plot(squeeze(stringLengthHistoryVert(1,2,:)))
plot(squeeze(stringLengthHistoryVert(2,2,:)))
plot(squeeze(stringLengthHistoryVert(3,2,:)))
plot(squeeze(stringLengthHistoryVert(4,2,:)))
grid on
title('Vertical Cable 2')
subplot(2,2,3)
hold on
plot(squeeze(stringLengthHistoryVert(1,3,:)))
plot(squeeze(stringLengthHistoryVert(2,3,:)))
plot(squeeze(stringLengthHistoryVert(3,3,:)))
plot(squeeze(stringLengthHistoryVert(4,3,:)))
grid on
title('Vertical Cable 3')
subplot(2,2,4)
hold on
plot(squeeze(stringLengthHistoryVert(1,4,:)))
plot(squeeze(stringLengthHistoryVert(2,4,:)))
plot(squeeze(stringLengthHistoryVert(3,4,:)))
plot(squeeze(stringLengthHistoryVert(4,4,:)))
grid on
title('Vertical Cable 4')
