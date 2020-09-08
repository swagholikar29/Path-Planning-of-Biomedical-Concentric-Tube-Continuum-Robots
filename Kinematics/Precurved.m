classdef Precurved < Robot
    % PRECURVED encapsulates the design parameters for a precurved tube
    % Original methods obtained from Endoscope class in https://github.com/comet-lab/Hale-Bopp on
    %   9/8/20020
    %   Original Author: Jesse F d'Almeida <jfdalmeida@wpi.edu>
    
    properties
        % Physical Endoscope
        OD                  % (m) outer diameter of endoscope
        precurve            % (1/m) precurved section curvature
        Ls                  % (m) length of straight section
        Lc                  % (m) length of curved section
        
        % Arc Parameters for forward kinematics
        kappa               % (1/m) curvature
        theta               % (rad) base rotation
        s                   % (m) arc length
        dz                  % (m) advancement
        
        % Forward Kinematics Parameters
        pose                % The position and orientation
        transformations     % Transformation matrix
        robotModel          % A model of the robot
        baseTransform       % base transformation (if any)
    end
    
    methods
        function self = Precurved(OD, precurve, Ls, Lc)
            % Class constructor
            %   only parameters that should change are current orientation
            
            self@Robot(2);      % parent class constructor, only 2 links: straight section and curved section
        
            self.OD = OD;
            self.precurve = precurve;
            self.Ls = Ls;
            self.Lc = Lc;
        end
        
        function self = fwkine(self, q, baseTransform)
            % Maps joint variable to arc parameters.
            %  fwkine(q, baseTransform) sets attributes pose and transform
            %  For now, assumes knowledge of arc parameters, since current
            %  joint variables are unknown
            %
            %   q = [theta, dz, k]    arc parameters of configuration
            %   theta [rad] = base rotation
            %   dz [m] = advancement
            %   k [1/m] = curvature   *optional for deformation* 
            %  baseTransform [4x4 matrix] *optional* homogenous transformation matrix
            %   where the endoscope begins
            
            % default baseTransform of none
            if ~exist('baseTransform', 'var')
                baseTransform = eye(4);
            end
            
            % parse the input q
            base_rotation = q(1);
            dz = q(2);
            if length(q) == 2
                curvature = self.precurve;
            else
                curvature = q(3);
            end
            
            % configurations for each section:
            %(base translations, curved, straight tip)
            kappa_list = [0 curvature];
            s_list = [dz self.Lc];
            theta_list = [base_rotation 0];
            
            % configurations
            base = [kappa_list(1) s_list(1) theta_list(1)];     % initial straight translation
            curve = [kappa_list(2) s_list(2) theta_list(2)];    % curved section
            c = [base curve]
            
            self.kappa = kappa_list;
            self.s = s_list;
            self.theta = theta_list;
            
            % forward kinematics based on configuration
            [P,T] = fwkine@Robot(self, c, baseTransform);
            self.pose = P(1:3,:);
            self.transformations = T;
            
            self.baseTransform = baseTransform;
        end
        
        function robotModel = makePhysicalModel(self)
            % Creates model of endoscope for plotting
            %  assumes constant curvature
            
            ptsPerM = 1e3;  % number of points for the curve
            
            P = self.pose;
            T = self.transformations;
            
            kappa = self.kappa;
            radius = 1 ./ kappa;
            s = self.s;
            
            backbone = P(:,1);     % 3d points of center
            
            for ii = 1 : size(P, 2) - 1
                if kappa(ii) == 0 % straight sections
                    
                    % generate points along a straight line
%                     distance = norm(P(:,ii+1) - P(:,ii));
%                     nPts = round(distance * ptsPerM);
                    % for straight section, only need 2 points
                    nPts = round(ptsPerM*1.5 * s(ii));
                    
                    X = linspace(P(1,ii),P(1,ii+1), nPts);
                    Y = linspace(P(2,ii),P(2,ii+1), nPts);
                    Z = linspace(P(3,ii),P(3,ii+1), nPts);
                    
                    backbone = [backbone [X(2:end);Y(2:end);Z(2:end)]];
                    
                else % curved sections
                    
                    % generate points along an arc of constant curvature
                    % and of length s
                    bend_angle = s(ii)*kappa(ii);
                    theta = linspace(0, bend_angle, s(ii)*ptsPerM);
                    
                    % points of the curve 
                    pts = radius(ii) .* [(1-cos(theta));
                                     zeros(1, length(theta));
                                     sin(theta);
                                     ones(1, length(theta)) / radius(ii)];
                    
                    backbone = [backbone ...
                        applytransform(pts(1:3,2:end), T(:,:,ii))]; 
                end
            end
            
            radius_vec = self.OD/2 * ones(1,size(backbone,2));   % radius for points along curve
            [X,Y,Z] = gencyl(backbone, radius_vec, 2, 10);       % gen cylinder around curve points
            
            % remove any NaNs
            X(any(isnan(X),2),:) = [];
            Y(any(isnan(Y),2),:) = [];
            Z(any(isnan(Z),2),:) = [];
            
            % set parameters of Model
            robotModel.backbone = backbone;
            robotModel.surface.X = X;
            robotModel.surface.Y = Y;
            robotModel.surface.Z = Z;
            
            self.robotModel = robotModel;
        end
        
        function t = getTransformations(self)
            t = self.transformations;
        end
    end
end