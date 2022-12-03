classdef Precurved < Robot
    % PRECURVED encapsulates the design parameters for a precurved tube
    % Original methods obtained from Endoscope class in https://github.com/comet-lab/Hale-Bopp on
    %   9/8/2020
    %   Original Author: Jesse F d'Almeida <jfdalmeida@wpi.edu>
    
    properties
        % Physical
        OD                  % (m) outer diameter
        ID                  % (m) inner diameter
        precurve            % (1/m) precurved section curvature
        Ls                  % (m) length of straight section
        Lc                  % (m) length of curved section

        % Material Properties
%         E = 70e6            % (Pa) Youngs Modulus (from wiki avg of austenite)
%         E = 1.83014e8       % (Pa) Young's Modulus for Bridge Nylon (https://taulman3d.com/bridge-nylon.html)
%         E = 1e9          % (Pa)
        E
        I                   % cross sectional moment of inertia of tube
        Poisson = 0.35      % Poisson's ratio for Nitinol
        
        % Arc Parameters for forward kinematics
        kappa               % (1/m) curvature
        phi               % (rad) base rotation
        s                   % (m) arc length
        dz                  % (m) advancement
        
        % Forward Kinematics Parameters
        pose                % The position and orientation
        transformations     % Transformation matrix
        robotModel          % A model of the robot
        baseTransform       % base transformation (if any)
        
        % Others
        handle              % plot model handle
    end
    
    methods
        function self = Precurved(OD, ID, precurve, Ls, Lc, E)
            % Class constructor
            %   only parameters that should change are current orientation
            
            self@Robot(2);      % parent class constructor, only 2 links: straight section and curved section
        
            self.OD = OD;
            self.ID = ID;
            self.precurve = precurve;
            self.Ls = Ls;
            self.Lc = Lc;
            self.E = E;
            
            self.I = (pi/32)*(OD^4-ID^4);
        end
        
        function self = fwkine(self, arcs, baseTransform)
            % Formats arc parameters to make transformations
            %  fwkine(q, baseTransform) sets attributes pose and transform
            %
            %   arc = [k, theta, dz]  matric of arc parameters of configuration
            %       k (1/m) = curvature   
            %       phi (rad) = rotations
            %       s (m) = lengths of each section
            %  baseTransform [4x4 matrix] *optional* initial transformation   
            %
            % default baseTransform of none
            if ~exist('baseTransform', 'var')
                baseTransform = eye(4);
            end
            
            numSections = size(arcs,1);
            self.nLinks = numSections;
            
            % parse the input q
            kappa_list = arcs(:,1)';
            phi_list = arcs(:,2)';
            s_list = arcs(:,3)';
            
            % add parameters to a single array for the fwkin
            c = [];
            for i = 1:numSections
                newC = [kappa_list(i) s_list(i) phi_list(i)];
                c = [c newC];
            end
            
            % store values
            self.kappa = kappa_list;
            self.s = s_list;
            self.phi = phi_list;
            
            % forward kinematics based on configuration
            [P,T] = fwkine@Robot(self, c, baseTransform);
            self.pose = P(1:3,:);
            self.transformations = T;
            
            self.baseTransform = baseTransform;
        end
        
        function robotModel = makePhysicalModel(self)
            % Creates model of endoscope for plotting
            %  assumes constant curvature
            
            ptsPerM = 2e3;  % number of points for the curve
            
            P = self.pose;
            T = self.transformations;
            
            kappa = self.kappa;
            radius = 1 ./ kappa;
            s = self.s;
            phi = self.phi;
            
            backbone = P(:,1);     % 3d points of center
            
            for ii = 1 : size(P, 2) - 1
                if kappa(ii) == 0 % straight sections
                    
                    % generate points along a straight line
%                     distance = norm(P(:,ii+1) - P(:,ii));
%                     nPts = round(distance * ptsPerM);
                    % for straight section, only need 2 points
                    nPts = ceil(ptsPerM*1.5 * s(ii));
                    
                    X = linspace(P(1,ii),P(1,ii+1), nPts);
                    Y = linspace(P(2,ii),P(2,ii+1), nPts);
                    Z = linspace(P(3,ii),P(3,ii+1), nPts);
                    
                    backbone = [backbone [X(2:end);Y(2:end);Z(2:end)]];
                    
                else % curved sections
                    
                    % generate points along an arc of constant curvature
                    % and of length s
                    bend_angle = s(ii)*kappa(ii);
                    nPts = ceil(ptsPerM * s(ii));
                    arcAng = linspace(0, bend_angle, nPts);
                    
                    % points of the curve 
                    pts = radius(ii) .* [(1-cos(arcAng)).*cos(phi(ii));
                                     (1-cos(arcAng)).*sin(phi(ii));
                                     sin(arcAng);
                                     ones(1, length(arcAng)) / radius(ii)];
                    
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