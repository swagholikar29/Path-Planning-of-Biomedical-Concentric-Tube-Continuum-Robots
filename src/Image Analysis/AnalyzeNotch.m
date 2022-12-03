function theta = AnalayzeNotch(notchImage,varargin)
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

theta = 0;
line_vec = zeros(2,2);
I = imshow(notchImage,'Parent',ax);
title(ax, "Select points on edge of notch cut");
for i = 1:2
    while(1)
        switch style
            case 'line'
                line = drawline('Color','magenta','Parent',ax);
                pos = line.Position;
            case 'points'
                point1 = drawpoint('Color','magenta','Parent',ax);
                point2 = drawpoint('Color','red','Parent',ax);
                pos = [point1.Position(1) point1.Position(2); 
                    point2.Position(1) point2.Position(2)];
                line = drawline('Position',pos,'Color','magenta');
                delete(point1); delete(point2);
        end
        line_vec(:,i) = [(pos(1,1) - pos(2,1));(pos(1,2) - pos(2,2))];
        line_vec(:,i) = line_vec(:,i)./norm(line_vec(:,i));
        choice = listdlg('PromptString',{'Are you happy with your line'},...
            'ListString',{'Yes','No'});
        if choice==1
            break;
        end
        delete(line)
    end
end
pause(0.1);

theta = acosd(dot(line_vec(:,1),line_vec(:,2))/(norm(line_vec(:,1))*norm(line_vec(:,2))));
end

