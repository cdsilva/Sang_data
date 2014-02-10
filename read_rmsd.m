function [t, rmsd] = read_rmsd(filename)

%%
fid = fopen(filename);

for i=1:5
    textscan(fid, '%s', 1, 'delimiter', '\n');
end

C = textscan(fid, '%f');

t = C{1}(1:2:end);
rmsd = C{1}(2:2:end);

fclose(fid);

