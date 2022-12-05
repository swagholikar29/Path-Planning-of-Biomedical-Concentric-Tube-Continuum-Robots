function qList = new_rrt(nPoints, delta_q, qBounds)
% qList = new_rrt(nPoints, delta_q, qBounds) performs rrt algorithm
%   nPoints = number of vertices
%   delta_q = a vector of the incremental distances for each parameter
%   qBounds = (opt.) a matrix [minBounds; maxBounds] of parameter bounds
%
%   qList = list of configurations
%       - normalized if qBounds is NOT provided
%       - scaled to bounds if qBounds is provided

close all
dim = length(delta_q);
qList = zeros(dim, nPoints);
hw = waitbar(0, 'Sampling the configuration space. Please wait...');

for n = 1 : nPoints - 1 
    q_rand = rand(dim,1); 
    q_near = nearestVertex(q_rand, qList, n);
    q_new = move(q_near, q_rand, delta_q);
    
    % Check if we are within bounds
    temp_delta_q = delta_q;
    while any(q_new  < 0) || any(q_new > 1)
        temp_delta_q = 0.5 * temp_delta_q; 
        q_new = move(q_near, q_rand, temp_delta_q);
    end
    
    qList(:,n+1) = q_new;
    waitbar(n/nPoints, hw, 'Sampling the configuration space. Please wait...');
end

% Scale if Given qBounds
if nargin == 3
    minBounds = qBounds(1,:)';
    maxBounds = qBounds(2,:)';
    scaleByBounds = @(q) q.*(maxBounds - minBounds) + minBounds;
    for index = 1:nPoints
        qList(:,index) = scaleByBounds(qList(:,index));
    end
end

close(hw);

%% Scatterplots and Histograms
% figure
% subplot(1, 3, 1)
% histogram(qListNormalized(1,:), 10) 
% title('Dimension 1')
% subplot(1, 3, 2)
% histogram(qListNormalized(2,:), 10)
% title('Dimension 2')
% if dim > 2
%     subplot(1, 3, 3)
%     histogram(qListNormalized(3,:), 10)
%     title('Dimension 3')
% end
% 
% if dim == 2
%     figure
%     scatter(qListNormalized(1,:), qListNormalized(2,:), 'b.')
%     title('Scatterplot of RRT Points')
%     xlabel('Dimension 1')
%     ylabel('Dimension 2')
%     axis equal
% else
%     figure
%     scatter3(qListNormalized(1,:), qListNormalized(2,:), qListNormalized(3,:), 'r.')
%     title('Scatterplot of RRT Points')
%     xlabel('Dimension 1')
%     ylabel('Dimension 2')
%     zlabel('Dimension 3')
%     axis equal 
end