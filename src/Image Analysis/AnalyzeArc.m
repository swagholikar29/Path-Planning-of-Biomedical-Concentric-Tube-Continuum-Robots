function radius_vec = AnalayzeArc(zoomedImg,varargin)
%ANALAYZENOTCH Has the user select 4 points to determine the angle the two
%cut sections are in relation to one another.
%
%   'Axis' - Optional Argument which is the axis to display the image one
%   'Style' - Name-Argument {'line','points'} which denotes if you want to
%   analyze a notch using lines or points.

%****** INPUT PARSING *********************
% default values
style = 'line';
styleOptions = {'line','points'};
numOfArcs = 1;

p = inputParser();
addRequired(p,'Image');
addOptional(p,'axis',0);
addParameter(p,'Style',@(x) any(validatestring(x,styleOptions)));
parse(p,path,varargin{:});

ax = p.Results.axis;
if ax == 0
    ax = gca;
end
style = p.Results.Style;
%*********************************************

radius_vec = [0 0];
line_vec = zeros(2,2);
I = imshow(zoomedImg,'Parent',ax);

for i = 1:1
    while(1)
        title(ax, "Select 3 points on the edge of the tube");
        switch style
            case 'line'
                line = drawpolyline('Color','magenta','Parent',ax);
                pos = line.Position;
            case 'points'
                point1 = drawpoint('Color','magenta','Parent',ax);
                point2 = drawpoint('Color','red','Parent',ax);
                point3 = drawpoint('Color','blue','Parent',ax);
                pos = [point1.Position(1) point1.Position(2);
                    point2.Position(1) point2.Position(2);
                    point3.Position(1) point3.Position(2)];
                line = drawpolyline('Position',pos,'Color','magenta');
                delete(point1); delete(point2); delete(point3);
        end
        
        choice = listdlg('PromptString',{'Are you happy with your curve'},...
            'ListString',{'Yes','No'});
        
        if choice ~=1
            delete(line);
            continue;
        end
        
        [radius cen] = fit_circle_through_3_points(pos);
        
        circle = drawcircle('Center', cen', 'Radius', radius);
        
        choice = listdlg('PromptString',{'Are you happy with your circle fit'},...
            'ListString',{'Yes','No'});
        
        radius_vec(i) = radius;
        delete(line)
        
        if choice==1
            break;
        else
            delete(circle);
        end
    end
end

pause(0.1);
end

