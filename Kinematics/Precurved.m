classdef Precurved < Robot
    % PRECURVED encapsulates the design parameters for a precurved tube
    % Originally obtained from https://github.com/comet-lab/Hale-Bopp on
    %   9/8/20020
    %   Original Author: Jesse F d'Almeida <jfdalmeida@wpi.edu>
    
    properties
        % Physical Endoscope
        OD = 5.38e-3         % [m] outer diameter of endoscope
        dcam = 1.25e-3       % [m] distance from center to focal point of camera
        dwrist = 1.25e-3     % [m] distance from center to center of wrist
        focal_height = 0    % [m] offset of focal point of camera from surface
        fov = deg2rad(85)   % [rad] camera field of view
        cam_range = 50e-3   % [m] range of camera view. Actual range is 3-50mm
        rbend_max = 25e-3   % [m] maximum bending radius
        bend_sec = 25.71e-3  % [m] length of bending section of endoscope
        tip_len = 11.43e-3   % [m] length of straight tip section
        
        % Arc Parameters for forward kinematics
        kappa               % [1/m] curvatute
        theta               % [rad] base rotation
        s                   % [m] arc length
        dz                  % [m] advancement
        
        % Forward Kinematics Parameters
        pose                % The position and orientation
        transformations     % Transformation matrix
        robotModel          % A model of the robot
        baseTransform       % base transformation (if any)
        
        % Other Transforms
        camT                % Transformation matrix from base to camera
        wristT              % Transformation matrix from base to wrist base
        
    end
    
    methods
        function self = Endoscope()
            % Class constructor
            %   only parameters that should change are current orientation
            
            self@Robot(3);      % parent class constructor, only 1 'link'
        end
        
        function self = fwkine(self, q, baseTransform)
            % Maps joint variable to arc parameters.
            %  fwkine(q, baseTransform) sets attributes pose and transform
            %  For now, assumes knowledge of arc parameters, since current
            %  joint variables are unknown
            %
            %   q = [k, theta, dz]    arc parameters of configuration
            %   k [1/m] = curvature
            %   theta [rad] = base rotation
            %   dz [m] = advancement
            %  baseTransform [4x4 matrix] *optional* homogenous transformation matrix
            %   where the endoscope begins
            
            % default baseTransform of none
            if ~exist('baseTransform', 'var')
                baseTransform = eye(4);
            end
            
            % parse the input q
            curvature = q(1);
            base_rotation = q(2);
            dz = q(3);
            
            % configurations for each section:
            %(base translations, curved, straight tip)
            kappa_list = [0 curvature 0];
            s_list = [dz self.bend_sec self.tip_len];
            theta_list = [base_rotation 0 0];
            
            % configurations
            base = [kappa_list(1) s_list(1) theta_list(1)];     % initial straight translation
            curve = [kappa_list(2) s_list(2) theta_list(2)];    % curved section
            tip = [kappa_list(3) s_list(3) theta_list(3)];      % straight tip section
            c = [base curve tip];
            
            self.kappa = kappa_list;
            self.s = s_list;
            self.theta = theta_list;
            
            % forward kinematics based on configuration
            [P,T] = fwkine@Robot(self, c, baseTransform);
            self.pose = P(1:3,:);
            self.transformations = T;
            
            center2cam = eye(4);             % transform from center of endoscope to camera
            center2cam(2,4) = -self.dcam;    % translated in neg y by cam offset
            self.camT = T(:,:,end)*center2cam;
            
            center2wrist = eye(4);              % transform from center of endoscope to wrist
            center2wrist(2,4) = self.dwrist;    % translated in pos y by wrist offset
            self.wristT = T(:,:,end)*center2wrist;
            
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
            
            backbone = P(:,1);     % 3d points of center of endoscope
            
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
            
            
            % == camera FoV cone ==
            n = 10;         % steps around cone
            h_step = linspace(0, self.cam_range, n);    % height
            th_step = linspace(0, 2*pi, n);             % theta around
            angle = tan(self.fov/2);                    % angle of triangle of vision
            
            [H,T] = meshgrid(h_step, th_step);
            
            % untranslated coords
            X = angle.*H.*cos(T);
            Y = angle.*H.*sin(T);
            Z = H;
            
            % initalize translated matrix
            Xtrans = zeros(n,n);
            Ytrans = zeros(n,n);
            Ztrans = zeros(n,n);
            
            % transform cone coords camera
            for i=1:n
                for j=1:n
                    point = [X(i,j); Y(i,j); Z(i,j)];
                    p = applytransform(point, self.camT);
                    Xtrans(i,j) = p(1);
                    Ytrans(i,j) = p(2);
                    Ztrans(i,j) = p(3);
                end
            end
            
            % add parameters to model
            robotModel.cam.X = Xtrans;
            robotModel.cam.Y = Ytrans;
            robotModel.cam.Z = Ztrans;
            
            self.robotModel = robotModel;
        end
        
        function inBounds = tipInBounds(self, Ptip)
            % Check if the wrist tip is out of bounds of the camera view
            % Ptip = [x,y,z] (mm) tip of wrist wrt base frame
            
            % transform tip position to camera frame
            ptip_cam = applytransform(Ptip, inv(self.camT));
            
            % check if tip is out of focal range
            if ptip_cam(3) > self.cam_range
                inBounds = false;
                return
            end
            
            u_ptip = ptip_cam / norm(ptip_cam); % unit vector of tip
            u_fov = [0 0 1]';   % unit vector for camera wrt camera
            
            % dot product of unit vector returns cos(angle between). Use
            % acos to get the angle between
            tip_angle = acos(dot(u_ptip, u_fov));
            
            % check if angle between is within fov angle
            if tip_angle > self.fov/2
                inBounds = false;
                return
            end
            
            inBounds = true;
        end
        
        function t = getTransformations(self)
            t = self.transformations;
        end
    end
end