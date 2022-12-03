function q_near = nearestVertex(q, qList, listLength)
    distances = pdist2(q', qList(:,1:listLength)');
    [~,i] = min(distances);
    
    if length(i) > 1
        disp('Multiple closest points detected');
    end
    
    q_near = qList(:,i);
end