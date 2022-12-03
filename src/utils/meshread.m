function mesh = meshread(filename)

fid = fopen(filename, 'rb');
mesh.id = fread(fid,1,'int32');
mesh.numverts = fread(fid,1,'int32');
mesh.numtris = fread(fid,1,'int32');
n = fread(fid,1,'int32');
if (n==-1)
    mesh.orient = fread(fid,3,'int32');
    mesh.dim = fread(fid,3,'int32');
    mesh.sz = fread(fid,3,'float');
    mesh.color = fread(fid,3,'int32');
else
    mesh.color = zeros(3,1);
    mesh.color(1) = n;
    mesh.color(2:3) = fread(fid,2,'int32');
end
mesh.vertices = reshape(fread(fid,3*mesh.numverts,'float'),[3,mesh.numverts]);
mesh.triangles = reshape(fread(fid,3*mesh.numtris,'int32'),[3,mesh.numtris]);
fclose(fid);