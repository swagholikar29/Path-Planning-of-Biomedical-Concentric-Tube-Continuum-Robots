%% Calculate tube diameters based on bending stiffness ratio
close all, clc, clear all;

%% anon functions
% second moment of inertia of tube based on outer and inner diameters
I = @(od, id) pi/32 * (od.^4 - (id).^4);

% ratio of bending stiffnesses of 2 tubes
clearance = .4;
getInnerID = @(d,w) (d-2*w-clearance);
% getInnerID = @(d,w) (0);
Mratio = @(d, w) I(d, d-w) / I(d-w-clearance, getInnerID(d,w));

% check if inner tube has negative inner diameter
isIDneg = @(d,w) getInnerID(d,w) < 0;

% calcs outer diameter of inner tube
getinnerOD = @(d,w) (d-w-clearance);


%% go through loops and calculate sets of tubes that meet specified ratio

wall = 1.2:.2:1.6;    % wall thicknesses
ods = 2:.2:6; % outer diameters 

% stores values of tubes corresponding to different ratios
%    under ratio of X - [outerOD  innerOD  wall  ratio]
Runder2 = [];
Rbetween2_3 = [];
Rbetween3_5 = [];
Rbetween5_8 = [];
Rover8  = [];

% loop through all possible combinations of tubes
for w = wall
    for d = ods
        % check if innermost diameter is negative
        if isIDneg(d,w), continue, end
        
        r = Mratio(d, w);   % ratio for this tube
        info = [d, getinnerOD(d,w), w, r];
        
        % sort tubes into respective matrix
        if r < 2
            Runder2(end+1,:) = info;
        elseif r < 3
            Rbetween2_3(end+1,:) = info;
        elseif r < 5
            Rbetween3_5(end+1,:) = info;
        elseif r < 8
            Rbetween5_8(end+1,:) = info;
        else
            Rover8(end+1,:) = info;
        end
    end
end
disp("Tubes under ratio of 2")
disp(Runder2)

disp("Tubes between ratio of 2 and 3")
disp(Rbetween2_3)

disp("Tubes between ratio of 3 and 5")
disp(Rbetween3_5)


disp("Tubes between ratio of 5 and 8")
disp(Rbetween5_8)





%% calculate diameters based on ratio
% w = 1.2;    % wall thickness diameter
% clearance = .4; % minimum clearance between tubes
% ratio = 4; % ratio of bending stiffnesses
% 
% syms d1 d2
% I1 = I(d1, d1-w);
% I2 = I(d2, d2-w);
% 
% I2 = subs(I2, d2, d1 - w - clearance);
% 
% eqn = I1 == ratio * I2;
% Ds = vpa(solve(eqn, 'MaxDegree', 4), 4);
% 
% fprintf("OD: %.2fmm    ID: %.2fmm\n", Ds(1), Ds(1)-w-clearance);




