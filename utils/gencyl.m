function [x y z] = gencyl( P, R , NBetween , NArround )
% gencyl
%
%   Generates a cylinder around a given centerline. 
%   Each point of the centerline must have an associated radius.
%
%   P is an array with the 3D centerline points column by column (3xn matrix).
%   R is a vector with the radii associated to each point of the centerline.
%   R should have as many elements as the columns of P.
%
%   NBetween is the number of 'steps' there should be between two centerline points.
%   It has to be >=2 (at least two points to join), and if it is >2, the radius values 
%   are automatically interpolated. So just give the centerline points and the exact 
%   same number of radius, and don't worry about the rest.
%
%   Narround is the number of equally spaced points arround the circumference of the 
%   cylinder. See 'doc cylinder'.
%   
%
%   Copyright 2011 Jonathan Hadida
%   ETS Montreal, Canada
%   Interventional Medical Imaging Laboratory
%
%   Contact: ariel dot hadida [a] google mail
    
    tovec = @(x)( x(:)' );
    N     = length(R);
    
    controlargs(nargin);
    
    % Directions
    D = diff( P , 1 , 2 );
    
    % Length between points
    L = sqrt(sum(D.^2));
    
    % Interpolate radius 
    if NBetween > 2
        R = interp1( 1:N , R , linspace(1,N, NBetween*(N-1) - (N-2) ) , 'spline' );
    end
    
    % Axis and angle of rotation for each cylinder
    % Careful, the order of the cross product is important for the rest
    Qvec = cross( [1e-6;1e-6;1]*ones(1,N-1) , D );
    Qang = acos(dot( [1e-6;1e-6;1]*ones(1,N-1) , D ) ./ L);
    
    % Allocation
    x = zeros( (N-1)*NBetween , NArround+1 );
    y = zeros( (N-1)*NBetween , NArround+1 );
    z = zeros( (N-1)*NBetween , NArround+1 );
    
    % Initialize loop indexes
    eidx = 0;
    erad = 1;
    
    for i=1:N-1
        
        % Start and end indexes
        sidx = eidx + 1;
        eidx = i*NBetween;
        
        % Radius start and end
        srad = erad;
        erad = srad + (NBetween-1);
        
        % Create cylinder
        [tx,ty,tz] = cylinder( R( srad:erad ) , NArround );
        
        % Rotation matrix
        Q = rotation_quat( Qvec(:,i) , Qang(i) );
        
        % Rotate and scale cylinder
        C = bsxfun( @plus , Q * vertcat( tovec(tx) , tovec(ty) , L(i)*tovec(tz) ) , P(:,i) );
        
        % Assign
        x( sidx:eidx , : ) = reshape( C(1,:)' , NBetween , NArround+1 );
        y( sidx:eidx , : ) = reshape( C(2,:)' , NBetween , NArround+1 );
        z( sidx:eidx , : ) = reshape( C(3,:)' , NBetween , NArround+1 );
        
    end
    
    function q = rotation_quat(vec,ang)
    % Quaternion matrix rotation
    
        cphi = cos(ang/2);sphi = sin(ang/2);
        vec  = vec / norm(vec);

        vals      = num2cell([ cphi tovec(vec)*sphi ]);
        [a b c d] = vals{:};

        q = [ 1 - 2*c^2 - 2*d^2 , 2*b*c - 2*d*a     , 2*b*d + 2*c*a ; ...
              2*b*c + 2*d*a     , 1 - 2*b^2 - 2*d^2 , 2*c*d - 2*b*a ; ...
              2*b*d - 2*c*a     , 2*c*d + 2*b*a     , 1 - 2*c^2 - 2*b^2 ];
    
    end

    function controlargs(argc)
        
        % 3D coordinates in col
        if size(P,1) ~= 3
            error('P should be a 3xn matrix.');
        end
        
        % Same number of points and radius
        if size(P,2) ~= N
            error('There should be the same number of points than radius.');
        end
        
        % At least two points
        if N < 2
            error('There should be at least 2 points.');
        end
        
        % Radii should be > 0
        if ~all(R > 0)
            error('Radius should be strictly positive scalars.');
        end
        
        % NBetween >= 1
        if argc < 3 || isempty(NBetween) || NBetween < 2
            NBetween = 2;
        else
            NBetween = ceil(NBetween);
        end
        
        % NArround > 2
        if argc < 4 || NArround < 3
            NArround = 20;
        else
            NArround = ceil(NArround);
        end
        
    end
    
end