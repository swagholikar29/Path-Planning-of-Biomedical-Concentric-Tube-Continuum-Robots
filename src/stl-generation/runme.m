% This script loads a stack of images of the vocal folds and converts them
% into an STL file

% Authors:        Emily Minch <evminch@wpi.edu>
%                 Rositsa Mihaleva <ramihaleva@wpi.edu>
%                 Alex Chiluisa <ajchiluisa@wpi.edu>
%                 Loris Fichera <lfichera@wpi.edu>
%
% Latest version: 6/15/2021

% Clean the workspace and close all the figures
clear, close all, clc

% What model are we working with?
modelID = 'L8C1'; % model ID
p = '\\research.wpi.edu\ROBOTICS\comet\NIDCD-2020\Vocal-Folds-Models-Bailly\LARYNX8\LARYNX8\L8a_pag_concaténé'; % path where the images are stored
prefix = 'L8a_pag_concat'; % image prefix
binarization_threshold = 29000;

% Read the image size
i = imread(fullfile(p, [prefix '0000.tif']));
X = size(i,2);
Y = size(i,1);
Z = size(dir(fullfile(p, '*.tif')), 1);

% Define the region of interest
t = readtable('image-rois.csv', 'ReadRowNames', true);
roi = t({modelID},:);
xLB = roi.XLowerBound; xUB = roi.XUpperBound;
yLB = roi.YLowerBound; yUB = roi.YUpperBound;

% Define the sampling interval across the three axes
dX = roi.SamplingX; dY=roi.SamplingY; dZ=roi.SamplingZ;

% Initialize the matrix where we are going to store all the images
VocalFolds = zeros((yUB - yLB)/dY, ...
                   (xUB - xLB)/dX,...
                   (Z/dZ));

% Load the images!
idz = 1;

f = waitbar(0,'Processing Images...');

for k = 0:dZ:(Z-1)
    idx = 1;
    idy = 1;
    
  % Generate the image file name
  if k>9 && k<100
      filename = [prefix '00' num2str(k) '.tif'];
  elseif k>99 && k<1000
      filename = [prefix '0' num2str(k) '.tif'];
  elseif k>999 && k<10000
      filename = [prefix num2str(k) '.tif'];
  else
      filename = [prefix '000' num2str(k) '.tif'];
  end

  raw_image = imread(fullfile(p, filename));  % Read the image
  cropped_image = raw_image(yLB:yUB,xLB:xUB); % Crop the image
  downscaled_image = zeros((yUB-yLB)/dY,(xUB-xLB)/dX); % Initialize the downscaled image
 
  % Downsample the image
  for r = 1 : dX : (xUB-xLB)
      for q = 1 : dY : (yUB-yLB)
          downscaled_image(idx,idy) = cropped_image(q,r);
          idx=idx+1; 
      end
      
      idx = 1;
      idy = idy+1;
  end
  
  % Binarize the image
  [L, W] = size(downscaled_image);
  BW = zeros(L,W);  
  BW = downscaled_image > binarization_threshold;
              
  VocalFolds(:,:,idz) = BW; % Save the result in the VocalFolds variable
  idz = idz + 1;
  
  waitbar(k/(Z-1),f,'Processing Images...');
end

close(f)

make_STL_of_Array([modelID '.stl'],VocalFolds,...
                  roi.VoxelSizeX * dX,...
                  roi.VoxelSizeY * dY,...
                  roi.VoxelSizeZ * dZ);