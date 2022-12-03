%% PLOTSIMS plots the total area (mm^2) of a set of simulations against the tested L and R.
%  looks at the csv given. Must have an L and R column.
% 
% Author: Phillip Abell <pbabell@wpi.edu>
%
close all, clear, clc

%% User changes these
% Change these values to be for the CSV you want to read.
fileToRead = 'Simulation-Larynx1-numberNotches.csv';

maxR = 18;
minR = 3;
deltaR = 5;

maxL = 20;
minL = 5;
deltaL = 5;

%% Should not need to change these
% Read in the table and check the results.
if ~exist(fileToRead, 'file')
    error('File %s does not exist', fileToRead)
end
dataTable = readtable(fileToRead, 'PreserveVariableNames', true);

if ismember('L', dataTable.Properties.VariableNames) == 0
    error('Field L does not exist in file %s', fileToRead)
end

if ismember('R', dataTable.Properties.VariableNames) == 0
    error('Field L does not exist in file %s', fileToRead)
end
Z = ones(4,4);
[lAxis, rAxis] = meshgrid(minL:deltaL:maxL,minR:deltaR:maxR);
for idx1 = 1:numel(lAxis(1,:))
    for idx2 = 1:numel(rAxis(:,1))
       Z(idx1, idx2) = table2array(dataTable(dataTable.L==lAxis(idx1,idx2) & dataTable.R==rAxis(idx1,idx2), {'VisibleArea(mm^2)'}));
    end
end

%% Create the Plot
s = surf(lAxis,rAxis, Z);
s.DataTipTemplate.DataTipRows(1).Label = 'L:';
s.DataTipTemplate.DataTipRows(2).Label = 'R:';
s.DataTipTemplate.DataTipRows(3).Label = 'Visible Area:';
xlim([5 20]);
ylim([3 18]);
zlabel('Visible Area mm^2');
xlabel('L mm');
ylabel('R mm');
