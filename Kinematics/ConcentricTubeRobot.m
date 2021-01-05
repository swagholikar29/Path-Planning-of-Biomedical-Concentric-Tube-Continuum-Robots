classdef ConcentricTubeRobot
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
    end
    
    methods
        function self = ConcentricTubeRobot(OD, ID, k, Ls, Lc, E)
            n = length(OD);
            tubes = zeros(1, n);
            
            for i = 1:n
                tubes(i) = Precurved(OD(i), ID(i), k(i), Ls(i), Lc(i), E(i));
            end
            
            self.OD = OD;
            self.ID = ID;
            self.k = k;
            self.Ls = Ls;
            self.Lc = Lc;
            self.E = E;
            self.I = [tubes.I];
        end
        
        function self = fwkine(q)
            %Calculates the forward kinematics of the complete robot
            %   converts joint to arcs to send to each tube
            %INPUT
            %   q = [nTubes x [translation rotation]] joint variables
            
            numOverlaps = 2*numTubes;

            % Calculate overlap lengths
            % matrix defining curved/straight/skip for each tube at each overlap
%             isCurved = zeros(numOverlaps, numTubes);
            [s isCurved] = calcLinkLengths(tubes, q(:,1));
    
            
        end
        
        %% -----HELPERS----
        function [k, phi] = inplane_bending(theta, isCurved)
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
        secExists = isCurved;
        secExists(secExists == -1) = 0;
        
        M = self.E .* self.I;
        
        theta = reshape(theta, size(M));
        
        kx = sum(M .* self.k .* cos(theta) .* isCurved) / sum(M * secExists);
        ky = sum(M .* self.k .* sin(theta) .* isCurved) / sum(M * secExists);
        
        k = sqrt(kx.^2 + ky.^2);
        phi = atan2(ky, kx);
        end
        
    end
end

