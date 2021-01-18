function [newImage] = RotateImage(origImage, varargin)
%ROTATEIMAGE Takes in an image and allows prompts user to rotate in 90deg
%increments
%   'Axis' - Optional Argument which is the axis to display the image one


%****** INPUT PARSING *********************
previousRegions = [];


p = inputParser();
addOptional(p,'axis',0);
parse(p,origImage,varargin{:});

ax = p.Results.axis;
if ax == 0
    ax = gca;
end
%****************************************

I = imshow(origImage,'Parent',ax);
newImage = origImage;

while(1)  
    choice = listdlg('PromptString',{'Rotate Image by 90deg?'},...
        'ListString',{'Yes','No'});
    
    if choice==1
        newImage = imrotate(newImage, 90);
        I = imshow(newImage,'Parent',ax);
    else
        break;
    end
    
end

end

