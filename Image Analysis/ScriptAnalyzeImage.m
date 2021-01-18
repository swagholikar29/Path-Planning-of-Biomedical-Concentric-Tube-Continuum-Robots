%% IMAGE ANALYSIS
% insert selected path and script use functions to find curvatures for
% tubes in each image on path
% STEPS OF SCRIPT:
%   1) Run script at desired path with specified params
%   2) Set scale of image using diameter of tube by selecting area around
%      tip up tube and selecting 2 points to form diameter
%   3) Select area where tube is curved to zoom in on
%   4) Select three points along inside edge of curvature for circle fit
%   5) Repeat for each tube
close all, clear, clc;
%% args to set

% path to folder with images
path = "C:\Users\jfd42\Documents\ThesisTrials\trials\PA12_newset\Tube_5b_3a\";
od = 3.8; % outer diameter

%% required vargs
relative = false;

imgtype = 'curvature';

%% analyze images
AnalyzeFolder(path, 'isRelative', false, 'ImgType', 'curvature', 'OD', 5.4)