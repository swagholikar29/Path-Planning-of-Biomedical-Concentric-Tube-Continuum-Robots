function RemovePoints(simulationID)
%% This script remove the reachable points that are outside of the 3D larynx model
%  The algorithm take the positions of the reachable points that lay
%  ooutside the larynx and remove from the list of positions to shoot the
%  laser beam.
%
%  Authors: A. Chiluisa <ajchiluisa@wpi.edu>
%           R. Tougas <rmtougas@wpi.edu>
%           E. Minch  <evminch@wpi.edu>
%           R. Mihaleva <ramihaleva@wpi.edu>
%  
% Last Version: 11/19/2021

fprintf('* Removing Points outside of the larynx *\n')
% add dependencies
addpath('kinematics')
addpath('utils')
addpath('utils/stlTools')
addpath('path-planning')
addpath('../anatomical-models')


load([simulationID '.mat']);
file2 = 'tissue_closed.stl';
testp = pList' * 1000;
testa = aList';
Stl2 = fullfile('..', 'anatomical-models', modelID, file2);
[vertices,faces,Name2] = stlRead(Stl2);
in = intriangulation(vertices,faces,testp);

% Plot figure results

figure('Name', [simulationID, '-Points Removed'])
subplot(1,2,1)
h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
set(h,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
hold on;
axis equal
plot3(testp(:,1),testp(:,2),testp(:,3),'b.');
plot3(testp(in==0,1),testp(in==0,2),testp(in==0,3),'ro');
title('Original Points')


subplot(1,2,2)
h1 = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
set(h1,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
hold on;
axis equal
plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'b.');
%plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');
title('Removed Points')

pList = ([testp(in==1,1),testp(in==1,2),testp(in==1,3)]')/1000;
aList = ([testa(in==1,1),testa(in==1,2),testa(in==1,3)]');

% Remove the Transformation matrices of pList

index2delete = find(in==0);
if ~isempty(index2delete)
    
    TListLen = length(TList); 
    index2deleteLen = length(index2delete); 
    TemTList = zeros(4,4,TListLen-index2deleteLen);%output
    TemqList = zeros(6,TListLen-index2deleteLen);
    TemqListNormalized = zeros(6,TListLen-index2deleteLen);
    TemxList = zeros(3,TListLen-index2deleteLen);
    jumper = index2deleteLen;

    for i = TListLen: -1: 1
        if (jumper > 0)
            if(~(i == index2delete(jumper)))
                TemTList(:,:,i-jumper) = TList(:,:,i);
                TemqList(:,i-jumper) = qList(:,i);
                TemqListNormalized(:,i-jumper) = qListNormalized(:,i);
                TemxList(:,i-jumper) = xList(:,i);
            else
                jumper = jumper - 1;
            end
        else
            TemTList(:,:,i) = TList(:,:,i);
            TemqList(:,i) = qList(:,i);
            TemqListNormalized(:,i) = qListNormalized(:,i);
            TemxList(:,i) = xList(:,i);
        end
    end
    
    TList = TemTList;
    qList = TemqList;
    qListNormalized = TemqListNormalized;
    xList = TemxList;   
end

save([simulationID '.mat'], '-v7.3');