function [newImage rectPosition] = SelectNotch(origImage,varargin)
%SELECTNOTCH Takes in an image and allows the user to select a region. The
%output is a zoomed in picture of that region.
%   'PreviousRegions' - Optional Argument is a nx4 array indicating previous
%   rectangles that have been analyzed where a row consists of
%   [xmin,ymin,xdistance,ydistance]. The previous regions will be selected
%   on the image. Default value is [].
%   'Axis' - Optional Argument which is the axis to display the image one


%****** INPUT PARSING *********************
previousRegions = [];
txt = "";


p = inputParser();
addRequired(p,'origImage',@isnumeric);
checkmat = @(x) isnumeric(x) && (size(x,2) == 4 || size(x,2) == 0);
addOptional(p, 'previousRegions', previousRegions, checkmat);
addOptional(p,'axis',0);
addOptional(p,'title',@isstring);
parse(p,origImage,varargin{:});

previousRegions = p.Results.previousRegions;
ax = p.Results.axis;
if ax == 0
    ax = gca;
end
txt = p.Results.title;
%****************************************

I = imshow(origImage,'Parent',ax);
title(txt);

% Display previously selected regions
for i= 1:size(previousRegions,1)
    rectangle('Position',previousRegions(i,:),'EdgeColor','red','LineWidth',1.5,'Parent',ax)
end

while(1)    
    % Select new region
    roi = drawrectangle('Parent',ax);
    rectPosition = roi.Position;
    xmin = round(roi.Position(1));
    ymin = round(roi.Position(2));
    xmax = round(roi.Position(1) + roi.Position(3));
    ymax = round(roi.Position(2) + roi.Position(4));
    newImage = origImage(ymin:ymax, xmin:xmax, :);
    
    choice = listdlg('PromptString',{'Confirm zoom'},...
        'ListString',{'Yes','No'});
    if choice==1
        break;
    end
    delete(roi);
end
% imshow(newImage)
end

