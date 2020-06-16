ref = '/temp/ulbl_p/PROFILES/33133_3000_TimeEvolTest/New/n.dat';
raw = load(ref);
rref = raw(:, 1);

dname = '/temp/ulbl_p/PROFILES/33133_3000_TimeEvolTest/Old2/';
files = dir(dname);
files = files(3:end);

for k = 1:numel(files)
    if (strcmp(files(k).name, 'Da.dat'))
        continue;
    end
    
    raw = load([dname, files(k).name]);
    raw = [rref, interp1(raw(:,1), raw(:,2), rref, 'spline', 'extrap')]';
    
    fid = fopen([dname, files(k).name], 'w');
    fprintf(fid, '%.15e\t%.15e\n',raw);
    fclose(fid);
end