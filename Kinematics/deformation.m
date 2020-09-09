function q = deformation(ks, arcs, thetas, d, baseRot, isCurved)
    % Maps the results of elastic deformation to a standard
    % configuration
    %   N = numSections - 1
    % INPUTS
    %   ks   = [Nx1] (1/m) array of curvatures
    %   arcs = [Nx1] (m) array of section arc lengths
    %   thetas = [Nx1] (rad) array of rotations
    %   d    = (m) initial straight translation
    %   baseRot = (rad) initial base rotation
    %   isCurved: [N] optional array mapping sections that are curved
    %       1 = curved, 0 = straight, -1 = skip
    % OUTPUT
    %   q = [kappa theta s] arc parameters
    
    numCurves = length(ks);
    
    % optional isCurved defined
    if ~exist('isCurved', 'var')
        isCurved = ones(1, numCurves);
    end

    numCurves = length(ks);

    q = zeros(numCurves + 1, 3);    % output
    
    % initial straight section, no curve, only translation
    q(1,:) = [0 baseRot d];
    
    % following curved sections (multiple due to elastic deformations
    for i = 1:numCurves
        k = ks(i);
        s = arcs(i);
        
        if isCurved(i) == -1
            k = 0;
            s = 0;
        end
            
        q(i+1,:) = [k thetas(i) s];
    end
end
