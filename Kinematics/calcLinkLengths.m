function [arcLengths, isCurved] = calcLinkLengths(tubes, p)
%CALCLINKLENGTHS Given a set of tubes and translations, calculate the
%resulting arc lengths of the links and which tubes belong to the links
%
%INPUTS
%   tubes = [N] array of Precurved tube objects
%   p = [N] array of translation arc parameter of tubes
%
%OUTPUTS
%   arcLengths = [2N] array of arclengths of links
%   isCurved = [2NxN] array mapping sections that are curved [Link x Tube]
%       1 = curved, 0 = straight, -1 = end

numTubes =  length(tubes);
numLinks = 2*numTubes;

% format outputs
arcLengths = zeros(1, numLinks);
isCurved = zeros(numLinks, numTubes);

Ls = zeros(1,numTubes); % straight section lengths
Lc = zeros(1,numTubes); % curved section lengths
for ii = 1:numTubes
    Ls(ii) = tubes(ii).Ls;
    Lc(ii) = tubes(ii).Lc;
end

%% Find absolute distances for each tube
dists = zeros(1, numLinks);

% distance of each transition point from origin
for ii = 0: numTubes-1
    dists(2*ii + 1) = p(ii + 1);            % curved section statrt
    dists(2*ii + 2) = p(ii + 1) + Lc(ii + 1);   % tube ends
end

%% Sort distances to calculate links

% sort distances
[sortedDists] = sort(dists);

% consecutive distances subtracted to get link lengths
arcLengths(1) = sortedDists(1) - 0;
for link = 2:numLinks
    arcLengths(link) = sortedDists(link) - sortedDists(link-1);
end

%% Figure out if tubes are straight, curved, or ended for each link

% for each tube, go through arc lengths and compare dists to Lc
for ii = 1:numTubes
    for ll = 1:numLinks
        d = sortedDists(ll);
        
        % default tube ended
        curvedness = -1;
        
        % trans point is before start of curved section: tube is straight
        if d <= p(ii)
            curvedness = 0;
            
        % trans point is after straight section but before end of tube:
        % tube is curved
        elseif d <= (p(ii) + Lc(ii))
            curvedness = 1;
        end
        
        % add result for link in tube
        isCurved(ll, ii) = curvedness; 
    end
end
end

