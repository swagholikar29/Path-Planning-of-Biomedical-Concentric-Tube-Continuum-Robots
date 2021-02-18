classdef ConcentricTubeRobot < handle
    %CONCENTRICTUBEROBOT is a subclass containing the entire kinematics of
    %the robot
    
    properties
        nTubes  % number of component tubes
        tubes   % nx1 vector of Precurved objects
        
        OD      % nx1 vector containing the tubes' inner diameters
        ID      % nx1 vector containing the tubes' inner diameters
        k       % nx1 vector containing the pre-curvature of each tube [mm^-1]
        Lc      % nx1 vector containing the lenghts of the precurved section of each tube
        Ls      % nx1 vector containing the lengths of straight section of each tube
        E       % nx1 vector containing the Young moduli of each tube
        I       % nx1 vector containins the second moment of inertia for each tube
        
        q       % nx[translation rotation] current configuration
    end
    
    methods
        function self = ConcentricTubeRobot(OD, ID, k, Ls, Lc, E)
            n = length(OD);
            self.nTubes = n;
            
            for i = 1:n
                tubes(i) = Precurved(OD(i), ID(i), k(i), Ls(i), Lc(i), E);
            end
            self.tubes = tubes;
            self.OD = OD;
            self.ID = ID;
            self.k = k;
            self.Ls = Ls;
            self.Lc = Lc;
            self.E = E;
            self.I = [tubes.I];
        end
        
        function fwkine(self, q)
            %Calculates the forward kinematics of the complete robot
            %   converts joint to arcs to send to each tube
            %INPUT
            %   q = [nTubes x [translation rotation]] joint variables
            self.q = q;
            numOverlaps = 2*self.nTubes;
            
            % Calculate overlap lengths
            % matrix defining curved/straight/skip for each tube at each overlap
            % isCurved = zeros(numOverlaps, numTubes);
            [s, isCurved] = calcLinkLengths(self.tubes, q(:,1));
            
            % get emergent curvatures for each overlap section
            newK = zeros(numOverlaps, 1);  % init new ks for each overlapped section
            phi = zeros(numOverlaps, 1);
            
            for  link = 1:numOverlaps
%                 [newK(link), phi(link)] = inplane_bending(self.tubes,q(:,2), isCurved(link,:));
                [newK(link), phi(link)] = self.inplane_bending(q(:,2), isCurved(link,:));
            end
            
            arcs = self.links2arcs(numOverlaps, newK, phi, isCurved, s);
            
            for i = 1:self.nTubes
                self.tubes(i).fwkine(arcs(:,:,i));
            end
        end
        
        %% -----HELPERS----
        function [k_eq, phi_eq] = inplane_bending(self, theta, isCurved)
            %INPLANE_BENDING uses the elastic deformation formulas found in Webster2009
            %  to calculate the resulting curvature
            %  INPUTS
            %   thetas   = [n] array of base rotations
            %   isCurved = [n] optional array mapping if this section is curved
            %       1 = curved, 0 = straight, -1 = skip
            %  OUTPUTS
            %   k   =  curvature of link
            %   phi =  rotation of link
            
            % array that modifies isCurved for if the tube is there
            % if section is straight, still exists, but if -1 then does not
            % exist
            secExists = isCurved;
            secExists(secExists == 0) = 1;
            secExists(secExists == -1) = 0;
            
            % array that modifies isCurved for numerator
            addCurvature = isCurved;
            addCurvature(addCurvature == -1) = 0;
            
            M = self.E .* self.I;
            
            theta = reshape(theta, size(M));
            
            kx =  sum(M .* self.k .* cos(theta) .* addCurvature) / sum(M .* secExists);
            ky = sum(M .* self.k .* sin(theta) .* addCurvature) ./ sum(M .* secExists);
            
            if isnan(kx), kx = 0; end
            if isnan(ky), ky = 0; end
            
            k_eq = sqrt(kx.^2 + ky.^2);
            phi_eq = atan2(ky, kx);
        end
        
        function arcs = links2arcs(self, numOverlaps, newK, phi, isCurved, Ls)
            arcs = zeros(numOverlaps, 3, self.nTubes);
            for t = 1:self.nTubes
                
                curr_phi = 0;     % keeps track of current angle
                for link = 1:numOverlaps
                    k = newK(link);
                    abs_phi = phi(link);
                    s = Ls(link);
                    
                    % update relative theta
                    if curr_phi ~= abs_phi
                        rel_phi = abs_phi - curr_phi;
                        curr_phi = curr_phi + rel_phi;
                    else
                        rel_phi = 0;
                    end
                    
                    % if the section ends
                    if isCurved(link,t) == -1 || s == 0
                        k = 0;
                        s = 0;
                    end
                    
                    arcs(link,:,t) = [k rel_phi s];
                end
            end
        end
        
        function psi = calc_torsional_flex(self, alpha)
            options = optimoptions('fsolve', 'OptimalityTolerance', 1e-10);
            dU = @self.energy_func;
            
            initials = self.q(:,2);
            psi = fsolve(dU, initials, options);
            
        end
        
        function F = energy_func(self, x)
            psi1 = x(1);
            psi2 = x(2);
            
            l1 = self.q(1,1);
            l2 = self.q(2,1);
            
            alpha1 = self.q(1,2);
            alpha2 = self.q(2,2);
            
            % consts
            J = 2 * self.I;
            
            v = .30; % .30 - .55
            G = self.E/(2*(1+v));
            
            c1 = G(1) * J(1)/self.Ls(1);
            c2 = G(2) * J(2)/self.Ls(2);
            c3 = self.E(1)*self.I(1)*self.E(2)*self.I(2)*self.k(1)*self.k(2)/...
                        (self.E(1)*self.I(1)+self.E(2)*self.I(2));
            
            f1 = c3*l2*sin(psi1 - psi2) - (c1*(2*alpha1 - 2*psi1))./2;
            f2 = -(c2*(2*alpha2 - 2*psi2))./2 - c3*l2*sin(psi1 - psi2);

            F = [f1; f2];
        end
        
        function plotTubes(self)
            
            h = zeros(1,self.nTubes);
            colors = distinguishable_colors(self.nTubes);
            figure('Name', 'Precurved Tubes');
            hold on
            
            for i = 1:self.nTubes
                model = self.tubes(i).makePhysicalModel();
                
                % mesh model of the tube
                h(i) = surf(model.surface.X, model.surface.Y, model.surface.Z,...
                    'FaceColor', colors(i,:));
                
                % create a triad (coord frame) for the transformations
                trans = self.tubes(i).transformations;
                for jj = 1:size(trans,3)
                    triad('Matrix', trans(:,:,jj), 'scale', 5e-3);
                end
                
                axis('image');
                view([135 30]);
                grid on;
                
                camlight('headlight');
                material('dull');
                axis equal
                xlabel('X (m)');
                ylabel('Y (m)');
                zlabel('Z (m)');
                title('Concentric Precurved Tubes with Deformation');
            end
            
        end
    end
end

