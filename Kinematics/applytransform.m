function P = applytransform(P,T)
%APPLYTRANSFORM Applies the homogenous transformation T to a list of points P. P
%must be a 3xn matrix, with each column representing a point in the
%Cartesian space.

P = [P; ones(1, size(P, 2))];
P = cellfun(@(p) T * p, num2cell(P, 1), 'UniformOutput', false);
P = cell2mat(P);
P = P(1:3,:);
end
