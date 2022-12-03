function data = getSimData(simulationID, filename)
%GETSIMDATA Retrieves and formats all the data from a simulation file to
%assemble into a csv. Optionally will append to given csv
%   loads the simulation and assembles all of the necessary data
%
sim = load([simulationID '.mat']);
    
d = sim.results;
data = {d.modelID d.nPoints d.dq d.runtime d.visFaces d.percFaces d.visArea d.percArea d.simID};

if nargin > 1
    if isfile(filename)
        fid = fopen(filename, 'A');
        fprintf(fid, '%s,%d,%.2f,%.2f,%d,%.2f,%.2f,%.2f,%s\n', data{:});
        fclose(fid);
    else
        header = {'ModelID' 'nPoints' 'DeltaQ' 'RunTime(min)' 'VisibleFaces' 'PercTotalFaces' 'VisibleArea(mm^2)' 'PercTotalArea' 'SimulationID'};
        fid = fopen(filename, 'w');
        fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s\n', header{:});
        fprintf(fid, '%s,%d,%.2f,%.2f,%d,%.2f,%.2f,%.2f,%s\n', data{:});
        fclose(fid);
    end
end
end

