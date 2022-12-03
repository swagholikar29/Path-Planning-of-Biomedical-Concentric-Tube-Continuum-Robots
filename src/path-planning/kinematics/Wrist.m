classdef Wrist < Robot
    % WRIST This class encapsulates the design parameters of a
    % notched-tube wrist.
    %
    %   Authors: Loris Fichera <lfichera@wpi.edu>
    %            Floris Van Rossum <fjvanrossum@wpi.edu>
    %
    %   Latest revision: 07/21/2019
    properties
        ID        % [m] tube inner diameter
        OD        % [m] tube outer diameter
        nCutouts  % [int] total number of cutouts
        cutouts   % struct array - each struct contains the following fields:
        %           [u - spacing between i-1 and i-th notch]
        %           [h - height of i-th cutout             ]
        %           [w - width of i-th cutout              ]
        %           [alpha - orientation of i-th cutout    ]
        
        % Forward Kinematics Parameters
        pose                % The position and orientation
        transformations     % Transformation matrix
        curvature           % Curvature of the wrist links
        arcLength           % Arc length of the wrist links
        robotModel          % A model of the robot
        
        ybar
        theta_max
        deltal_max
    end
    
    methods
        function self = Wrist(ID, OD, nCutouts, cutouts)
            % Invoke the constructor of the parent class
            nLinks = nCutouts * 2 + 1;
            self@Robot(nLinks);
            
            % Copy geometric design parameters in local attributes
            self.ID = ID;
            self.OD = OD;
            self.nCutouts = nCutouts;
            self.cutouts = cutouts;
            
            % Pre-calculate some variables to make the forward kinematics run faster
            ro = self.OD / 2; % outer radius of tube in [m];
            ri = self.ID / 2; % inner radius of tube in [m];
            
            h = self.cutouts.h; % Height of the cutouts in [m]
            w = self.cutouts.w; % Cut depth in [m]. See Figure 4 again.
            d = w-ro; % intermediate variable. Depth of cut as measured from y = 0. See Figure 4.
            
            
            for ii = 1 : nCutouts
                phio = 2 * acos(d(ii) / ro);
                phii = 2 * acos(d(ii) / ri);
                
                ybaro = (4 * ro * (sin(0.5 * phio)) ^ 3)/ (3 * (phio - sin(phio)));
                ybari = (4 * ri * (sin(0.5 * phii)) ^ 3)/ (3 * (phii - sin(phii)));
                
                Ao = ( (ro ^ 2) * ( phio - sin(phio))) / 2;
                Ai = ( (ri ^ 2) * ( phii - sin(phii))) / 2;
                
                self.ybar(ii) = (ybaro * Ao - ybari * Ai) / (Ao - Ai);
                
                % Calculate the maximum bending for this cutout
                self.theta_max(ii) = h(ii) / (ro + self.ybar(ii));
                
                % Calculate the tendon displacement for this cutout when
                % this hits the hard stop.
                self.deltal_max(ii) = h(ii) - (ro - ri) * self.theta_max(ii);
            end
        end
        
        function cutoutDispl = calcdispl(self, totalDisplacement)
            % This function calculates the tendon displacement experienced
            % by each individual cutout. By default, the total tendon
            % displacement is simply uniformly redistributed among the
            % existing cutouts. When a cutout hits the hard stop, it is no
            % longer articulable.
            
            % First assign each cutout 1/n of the total displacement
            cutoutDispl = totalDisplacement ./ self.nCutouts .* ...
                ones(1, self.nCutouts);
            
            overHardStop = cutoutDispl > self.deltal_max;
            
            if any(overHardStop)
                while true
                    % calculate the excess displacement on the cutouts that have
                    % hit the full stop
                    excessDisplacement = sum((cutoutDispl - self.deltal_max) .* double(overHardStop));
                    
                    % redestribute among the other cutouts
                    cutoutDispl(overHardStop) = self.deltal_max(overHardStop);
                    cutoutDispl(~overHardStop) = cutoutDispl(~overHardStop) + ...
                        excessDisplacement ./ sum(~overHardStop);
                    
                    % check if any cutouts are still running over the hard
                    % stop
                    overHardStop = cutoutDispl > self.deltal_max;
                    
                    if all(cutoutDispl >= self.deltal_max)
%                         warning('All the cutouts hit the hard stop');
                        cutoutDispl = self.deltal_max;
                        break;
                    elseif excessDisplacement < 1e-9
                        break;
                    end
                end
            end
        end
        
        function c = jointvar2arcparams(self, q)
            % == Get the endoscope joint variables q
            tendonDisplacement = q(1); % tendon displacement [m]
            t_rot   = q(2);            % tube rotation [m]
            t_adv   = q(3);            % tube advancement [m]
            
            cutoutDispl = self.calcdispl(tendonDisplacement);
            
            % == Read the geometric design parameters
            ri = self.ID / 2; % inner radius of tube in [m];
            
            h = self.cutouts.h; % Height of the cutouts in [m]
            
            % == Perform a robot-specific mapping from the joint variables q
            % to the vector of arc parameters c
            c = zeros(1, 3 * self.nLinks);
   
            % initialize a counter to iterate on the notches
            kk = 1;
            
            self.curvature = zeros(1,self.nLinks);
            self.arcLength = zeros(1,self.nLinks);
            
            % keep track of length wrist upon reaching each notch
            length_at_section = 0;
            
            for jj = 0 : self.nLinks - 1
                % Even-numbered links are notches
                % Odd-numbered links are straight sections
                % Links lower than displacement are also straight sections
                if mod(jj+1,2) == 1
                    % For a straight section, the curvature is always zero
                    kjj = 0;
                    if jj+1 == 1 && t_adv >= 0
                        % For the first staight section only, the length is
                        % given by the robot advancement if positive, and
                        % the rotation is given by the tube rotation
                        sjj = t_adv;
                        thetajj = t_rot;
                    elseif jj+1 == 1 && t_adv < 0
                        % If the advancement is negative, the first section
                        % is pure rotation.
                        sjj = 0;
                        thetajj = t_rot;
                    else
                        % For all successive uncut sections, there is no
                        % rotation, and the advancement is given by the
                        % length of the uncut section itelf
                        sjj = self.cutouts.u(kk-1); % Uncut section [m]
                        thetajj = self.cutouts.alpha(kk-1);
                        if length_at_section < abs(t_adv) && t_adv < 0
                            % If the length of the wrist at this point is
                            % less than that of the advancement when it the
                            % advancement is negative, do not include for
                            % visualization.
                            if length_at_section + sjj >= abs(t_adv)
                                % if part of the section will appear
                                % outside of the endoscope, include that
                                % portion
                                og_length_at_section = length_at_section;
                                length_at_section = length_at_section + sjj;
                                sjj = sjj + og_length_at_section - abs(t_adv);
                            else
                                length_at_section = length_at_section + sjj;
                                sjj = 0;
                                kjj = 0;
                            end
                        end
                    end
                    
                else
                    % For a curved section (i.e. a notch), calculate the
                    % arc parameters using the relations in [Swaney2017]
                    kjj = (cutoutDispl(kk)) / (h(kk) * (ri + self.ybar(kk)) - cutoutDispl(kk) * self.ybar(kk));
                    sjj = h(kk) / ( 1 + self.ybar(kk) * kjj); % original kappa
                    thetajj = 0;
                    if length_at_section < abs(t_adv) && t_adv < 0
                        if length_at_section + sjj >= abs(t_adv)
                            og_length_at_section = length_at_section;
                            length_at_section = length_at_section + sjj;
                            sjj = sjj + og_length_at_section - abs(t_adv);
                        else
                            sjj = self.cutouts.h(kk);
                            length_at_section = length_at_section + sjj;
                            sjj = 0;
                            kjj = 0;
                        end
                    end
                    % move to the next cutout
                    kk = kk + 1;
                end
                
                c(jj*3+1:jj*3+3) = [kjj sjj thetajj];
                
                self.curvature(jj+1) = kjj;
                self.arcLength(jj+1) = sjj;
            end
        end
        
        function fwkine(self, q, baseTransform)
            % map joint variables to arc parameters
            c = self.jointvar2arcparams(q);
            
            % Calculate the robot-independent mapping and return the result
            [P,T] = fwkine@Robot(self, c, baseTransform);
            
            % Save pose and transformations in local attributes
            self.pose = P(1:3,:);
            self.transformations = T;
            self.transformations(1:3,4,:) = self.transformations(1:3,4,:);
        end
        
        
        function [Jrobot,Jp] = jacob0(self, q)
            % Read the joint variables
            t_displ = q(1); % tendon displacement [m]
            
            % Read the geometric design parameters
            ri = self.ID / 2; % inner radius of tube in [m];
            h = self.cutouts.h; % Height of the cutouts in [m]
            
            % Map joint variables to arc parameters
            c = self.jointvar2arcparams(q);
            
            % Calculate the robot-independent Jacobian
            J = jacob0@Robot(self, c);
            
            % Expand the time derivatives of the arc parameters using the
            % chain rule
            cdot = zeros(length(c), 1);
            cdotmatrix = zeros(3*self.nLinks, 3);
            
            % initialize a counter to iterate on the notches
            kk = 1;
            
            for jj = 0 : self.nLinks - 1
                % Even-numbered links are notches
                % Odd-numbered links are straight sections
                
                if mod(jj+1,2) == 1 % for a straight section
                    kdotjj = 0;
                    
                    if jj+1 == 1
                        % For the first staight section only, the rate of
                        % change of advancement and rotation are directly
                        % controlled through manipulation of the joint
                        % variables
                        sdotjj     = 1;
                        thetadotjj = 1;
                    else
                        % For all successive uncut sections, there is no
                        % translational nor rotational velocity;
                        sdotjj     = 0;
                        thetadotjj = 0;
                    end
                else % for a notched section
                    %                     d = h(kk) * (ri + self.ybar(kk)) - t_displ * self.ybar(kk);
                    %                     kdotjj = 1/d + self.ybar(kk) *  t_displ / d^2;
                    %                     sdotjj = -(h(kk) * self.ybar(kk) / d + self.ybar(kk)^2 * t_displ / d^2) * ...
                    %                         1 / (1 + (t_displ * self.ybar(kk)) / d)^2;
                    kdotjj = 1/(h(kk)*(ri + self.ybar(kk)) - self.ybar(kk)*t_displ) + (self.ybar(kk)*t_displ)/(h(kk)*(ri + self.ybar(kk)) - self.ybar(kk)*t_displ)^2;
                    sdotjj = -(h(kk)*((self.ybar(kk))/(h(kk)*(ri + self.ybar(kk)) - self.ybar(kk)*t_displ) + (self.ybar(kk)^2*t_displ)/(h(kk)*(ri + self.ybar(kk)) - self.ybar(kk)*t_displ)^2))/((self.ybar(kk)*t_displ)/(h(kk)*(ri + self.ybar(kk)) - self.ybar(kk)*t_displ) + 1)^2;
                    thetadotjj = 0;
                    
                    % move on to the next cutout
                    kk = kk + 1;
                end
                
                cdot(jj*3+1:jj*3+3) = [kdotjj sdotjj thetadotjj];
                cdotmatrix(jj*3+1,1) = kdotjj;
                
                if jj+1 == 1
                    cdotmatrix(jj*3+2,3) = sdotjj;
                else
                    cdotmatrix(jj*3+2,1) = sdotjj;
                end
                
                cdotmatrix(jj*3+3,2) = thetadotjj;
            end
            
            % Calculate the forward kinematics - required to finalize
            % calculation of the Jacobian
            self.fwkine(q, eye(4));
            pCurrent = self.pose(:,end); % the forward kinematics returns values in mm
            
            Jrobot = J*cdotmatrix;
            Jp = Jrobot(1:3,:) - skew(pCurrent) * Jrobot(4:6,:);
        end
        
        
        function q = invkine(self, pTarget, qCurrent)
            pTarget = pTarget;
            
            while true
                % calculate the current location
                self.fwkine(qCurrent, eye(4));
                pCurrent = self.pose(1:3,end);
                
                % calculate the difference between current and target location
                err = norm(pTarget - pCurrent)
                
                % if the difference < epsilon, return
                if err < 1e-1, break; end
                
                % else, update the "joint" variables using the inverse of the
                % jacobian
                [~,Jp] = self.jacob0(qCurrent);
                %Jp = J(1:3,:) - skew(pCurrent) * J(4:6,:)
                Jpinv = pinv(Jp);
                
                K = diag([100 100 100]);
                deltaQ = K * Jpinv * (pTarget - pCurrent);
                deltaQ(1) = deltaQ(1);
                deltaQ(3) = deltaQ(3);
                qCurrent = qCurrent + deltaQ';
                
                RM = self.makePhysicalModel();
                h.XData = RM.surface.X;
                h.YData = RM.surface.Y;
                h.ZData = RM.surface.Z;
                
                drawnow
                
                %
                %                 triad('Matrix', eye(4), 'linewidth', 2.5);
                %                 triad('Matrix', self.transformations(:,:,end), 'linewidth', 2.5);
                %
                %                 RM = self.makePhysicalModel();
                %
                %                 X = RM.surface.X;
                %                 Y = RM.surface.Y;
                %                 Z = RM.surface.Z;
                %                 surf(X, Y, Z, 'FaceColor','green');
                %                 hold on
                %                 drawnow
            end
        end
        
        
        function robotModel = makePhysicalModel(self)
            ptsPerM = 3e3;
            
            P = self.pose;
            T = self.transformations;
            
            kappa = self.curvature;
            radius = 1 ./ kappa;
            s = self.arcLength;
            
            robotBackbone = P(:,1);
            
            for ii = 1 : size(P, 2) - 1
                if kappa(ii) == 0
                    
                    % generate points along a straight line
                    distance = norm(P(:,ii+1) - P(:,ii));
                    nPts = round(distance * ptsPerM);
                    
                    X = linspace(P(1,ii),P(1,ii+1), nPts);
                    Y = linspace(P(2,ii),P(2,ii+1), nPts);
                    Z = linspace(P(3,ii),P(3,ii+1), nPts);
                    
                    robotBackbone = [robotBackbone [X(2:end);Y(2:end);Z(2:end)]];
                    %                     height_at_section = height_at_section + self.cutouts.u(ii);
                    
                else % cutouts: curved sections & anything over the advancement
                    
                    % generate points along an arc of constant curvature
                    % and of length s
                    
                    % !FIXME ensure ptspermm is met
                    theta = 0 : s(ii)*kappa(ii)/5 : s(ii)*kappa(ii);
                    
                    pts = radius(ii) .* [(1-cos(theta));
                        zeros(1, length(theta));
                        sin(theta);
                        ones(1, length(theta)) / radius(ii)];
                    
                    robotBackbone = [robotBackbone ...
                        applytransform(pts(1:3,2:end), T(:,:,ii))];
                end
            end
            
            %robotBackbone = applytransform(robotBackbone, baseTransform);
            
            radiusVec = self.OD/2*ones(1,size(robotBackbone,2));
            [X,Y,Z] = gencyl(robotBackbone, radiusVec, 2, 10);
            
            robotModel.backbone = robotBackbone;
            robotModel.surface.X = X;
            robotModel.surface.Y = Y;
            robotModel.surface.Z = Z;
            self.robotModel = robotModel;
        end
    end
end