classdef Robot < handle
    % Robot This class implements the robot-independent kinematics of a continuum robot.
    %   Author: Loris Fichera <lfichera@wpi.edu>
    %   Latest revision: 07/07/2019

    properties
        nLinks % number of links
    end
    
    methods
        function self = Robot(nLinks)
            self.nLinks = nLinks;
        end
        
        function [P, T] = fwkine(self, c, baseTransform)
            % c is the configuration vector, i.e. the vector that contains
            % the arc parameters for each link in the following order:
            % curvature, arc length, base rotation (k, s, theta)
            %
            % baseTransform is an optional parameter that defines a
            % transformation to be applied to the base of the tube -
            % useful, for instance, to register the workspace of the robot
            % with another space
            
            % == Iteratively calculate the forward kinematics ==
            % 1. Initialize pose and transformation matrices
            T = repmat(eye(4), 1, 1, self.nLinks + 1);
            T(:,:,1) = baseTransform;
            
            P = zeros(4, self.nLinks + 1);
            P(1:3,1) = baseTransform(1:3,end);
            P(4,:) = ones(1, self.nLinks + 1);
            
            % 2. Using the product of exponentials formula, generate the
            % transformations for each link and propagate the
            % transformations
            for jj = 0 : self.nLinks - 1
                cjj = c(jj*3+1:jj*3+3);                
                T(:,:,jj+2) = T(:,:,jj+1) * arckinematics(cjj);
                P(1:3,jj+2) = T(1:3,end,jj+2);
            end
        end
        
        % this function calculates the Jacobian for a single link
        function J = jacobiansinglelink(~,k,s,theta)
            if k == 0 % when curvature is null, need to apply the limit for k -> 0
                J = [cos(theta) * -s^2 * 0.5  0 0;
                     sin(theta) * -s^2 * 0.5  0 0;
                     0                        1 0;
                    -s * sin(theta)           0 0;
                     s * cos(theta)           0 0;
                     0                        0 1];
            else
                J = [cos(theta) * (cos(k * s) - 1) / k^2 0              0;
                    sin(theta) * (cos(k * s) - 1) / k^2  0              0;
                    -(sin(k*s) - k * s) / k^2            1              0;
                    -s * sin(theta)                     -k * sin(theta) 0;
                    s * cos(theta)                       k * cos(theta) 0;
                    0                                    0              1];
            end
        end
        
        
        % this function calculates the robot-independent Jacobian
        function J = jacob0(self, c)               
             % initialize the Jacobian matrix
             J = zeros(6, 3 * self.nLinks);
             
             % Add the Jacobian for the first link
             J(:,1:3) = self.jacobiansinglelink(c(1), c(2), c(3));
             
             % iteratively calculate and stack the jacobians for the other
             % links
             T = eye(4);
             
             for ii = 0 : self.nLinks - 2
                 % calculate the Adjoint transformation first
                 Tii = arckinematics(c(ii*3+1:ii*3+3));
                 T = T * Tii;
                 Ad = [T(1:3,1:3), skew(T(1:3,end)) * T(1:3,1:3);
                         zeros(3), T(1:3,1:3)];
                     
                 % then calculate the jacobian for the current link
                 Jii = self.jacobiansinglelink(c(ii*3+4), c(ii*3+5), c(ii*3+6)); 
                 J(:,ii*3+4:ii*3+6) = Ad * Jii;
             end
        end
    end
end