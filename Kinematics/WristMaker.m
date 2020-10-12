classdef WristMaker < handle
    %WRISTMAKER Provides helper classes to analyze and create wrists
    %   
    %   Original Author: Jesse F d'Almeida
    %   Last revision: 10/11/2020
    
    properties
        % wrist and key parameters
        robot       % Wrist Object
        ybar        % (m) neutral bending plane
        OD          % (m) outer diameter
        ID          % (m) inner diameter
        ro          % (m) outer radius
        ri          % (m) inner radius
        yeild       % (Pa) yield strength of material
        modulus     % (Pa) Youngs modulus of material
        emax        % yield strain of material
        Ao          % Second moment of inertia of outer section
        Ai          % second moment of inertia of inner section
        
        % notch configuartions
        n           % number of notches
        g           % (m) depth of notch cut
        u           % (m) uncut section height
        h           % (m) cut section height
        
        % notch properties
        kj
        thetaj
        sj
        
        % calculated properites
        max_theta   % (rad) max bending radius of tube
        max_thetaj  % (rad) max bending radius of single notch
        s
        
        
    end
    
    methods
        function self = WristMaker(ID, OD, n, h, u, g)
            % assign properties
            self.OD = OD;
            self.ro = OD/2;
            self.ID = ID;
            self.ri = ID/2;
            self.n = n;
            self.h = h;
            self.u = u;
            self.g = g;
            
            % Create Wrist
            cutouts.w = g * ones(1,n); % [m]
            cutouts.u = u * ones(1,n); % [m]
            cutouts.h = h * ones(1,n); % [m]
            cutouts.alpha = zeros(1,n);
            
            self.robot = Wrist(ID, OD, n, cutouts);
            
            % get wrist params
            self.ybar = self.robot.ybar(1);
            self.Ao = self.robot.Ao(1);
            self.Ai = self.robot.Ai(1);
        end
        
        function [theta, R, s] = maxBendFromGeometry(self)
            % maxBendFromStain calculates the maximum possible bending
            % before permenant deformation based on the wrist material
            % OUTPUTS
            %   theta = (rad) maximum bending angle
            %   R     = (m)   minimum bending radius
            %   s     = (m)   arc length of curved section
            
            theta = self.n * self.h / (self.ro + self.ybar);
            self.max_theta = theta;
            
            R = self.ro + (self.n - 1) * self.u / theta;
            s = self.h * self.n + self.u * (self.n - 1);
        end
        
        function [theta, k, s] = maxBendFromStain(self)
            % maxBendFromStain calculates the maximum possible bending
            % before permenant deformation based on the wrist material
            % OUTPUTS
            %   theta = (rad) maximum bending angle
            %   k     = (1/m) maximum curvature
            %   s     = (m)   arc length of curved section
            y = [self.ro - self.ybar, self.ybar - (self.g - self.ro)];
            [d, idx] = max(abs(y));   
            d = d * sign(y(idx));

            kj = self.emax/(d - self.ybar * self.emax);    % calc max kappa

            thetaj = (self.h/(1 + self.ybar * kj))*kj;    % bending angle of each notch
            theta = thetaj * self.n;                 % total bending angle of tube

            sj = self.h / (1 + self.ybar * kj);    % arc length of each notch
            s = sj * self.n + self.u * (self.n-1);

            % calc min row
            rho = s/theta;
            k = 1/rho;
            
            self.kj = kj;
            self.thetaj = thetaj;
            self.sj = sj;
        end
        
        function [Emin, stress, strain] = calcNotchForce(self)
            Rj = 1/kj;
            
            Fmax = Rj * self.modulus * (self.Ao - self.Ai) / (self.ro + self.ybar)
    
            notch_area = pi * (self.ro^2 - self.ri^2) - (self.Ao - self.Ai);   % surface area of notch

            stress = .5 * Fmax * self.g / notch_area

            strain = (self.h - self.thetaj * (Rj + self.ro))/(self.thetaj * (Rj + self.ro));

            Emin = stress / strain;
            
        end
    end
end

