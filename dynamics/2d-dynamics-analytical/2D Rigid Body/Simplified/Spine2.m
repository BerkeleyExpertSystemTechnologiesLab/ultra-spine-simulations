classdef Spine2 < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lSize = [];
        Angle = [];
    end
    
    properties (SetAccess = private)
        bars = struct('lx',[],'ly',[],'theta',[]);
        wires = struct('rx',[],'ry',[],'thr',[]);%thr stands for theta r
        x1 = [0 0];
        x2 = [0 2];
        th = [];
    end
    
    properties (Hidden)
       k = 30;
       ro = 0.3;
    end
    
    methods
        function obj = Spine2(a,b)
            obj.lSize = a;
            obj.Angle = b;
            obj.buildStruct;
%             obj.plot;
        end
        
        function plot(obj)
           hold on;
           for n = 1:3

              l1x = obj.bars(n).lx;
              l1y = obj.bars(n).ly;
              
              l2x = obj.bars(n+3).lx;
              l2y = obj.bars(n+3).ly;
              


              plot(l1x,l1y,'b');
              plot(l2x,l2y,'r');
              
           end
           
           for n = 1:5
               r1x = obj.wires(n).rx;
               r1y = obj.wires(n).ry;
               plot(r1x,r1y,'k--');
           end
        end
    end
    
    methods (Access = private)
        function buildStruct(obj)
            lengths = obj.lSize;
            angle = obj.Angle;
            
            theta(1) = (180-angle)/2;
            theta(2) = theta(1) + angle;
            theta(3) = 270;
            obj.th = theta.*pi./180;
            X1 = obj.x1;
            X2 = obj.x2;
            
            for n = 1:3
%                 keyboard
                obj.bars(n).lx = [X1(1);
                    X1(1) + lengths*cosd(theta(n))];
                obj.bars(n).ly = [X1(2);
                    X1(2) + lengths*sind(theta(n))];
                obj.bars(n+3).lx = [X2(1);
                    (X2(1)+lengths*cosd(theta(n)))];
                obj.bars(n+3).ly = [X2(2);
                    (X2(2)+lengths*sind(theta(n)))];
                
                obj.wires(n).rx = [obj.bars(n).lx(2);
                    obj.bars(n+3).lx(2)];
                obj.wires(n).ry = [obj.bars(n).ly(2);
                    obj.bars(n+3).ly(2)];
                
%                 keyboard
                
                obj.bars(n).theta = theta(n);
                obj.bars(n+3).theta = theta(n);
            end
            
            obj.wires(4).rx = [obj.bars(2).lx(2);obj.bars(6).lx(2)];
            obj.wires(5).rx = [obj.bars(1).lx(2);obj.bars(6).lx(2)];
            
            obj.wires(4).ry = [obj.bars(2).ly(2);obj.bars(6).ly(2)];
            obj.wires(5).ry = [obj.bars(1).ly(2);obj.bars(6).ly(2)];
%             keyboard
            for i = 1:5
                dx = diff(obj.wires(i).rx);
                dy = diff(obj.wires(i).ry);
                obj.wires(i).thr = asind(dy/dx); 
            end
            
%             keyboard
        
        end
        
    end
    
end

