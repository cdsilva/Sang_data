function [pos, t, residue_ind, isH, isCA, box] = read_gro_file(filename)

%%
fid = fopen(filename);
nlines = 0;
while(~feof(fid))
    textscan(fid, '%s', 1, 'delimiter', '\n');
    nlines = nlines + 1;
end
fclose(fid);

%%

dim = 3;
fid = fopen(filename);

C = textscan(fid,'%s %f',1,'delimiter','=');
tmp_t = C{2};

C = textscan(fid, '%d', 1);
natoms = C{1};

nsteps = nlines / (natoms + 3);
if mod(nsteps, 1) ~= 0
    disp('ERROR');
    return;
end

t = zeros(nsteps, 1);
box = zeros(nsteps, dim);
pos = zeros(natoms, dim, nsteps);
residue_ind = zeros(natoms, 1);
isCA = false(natoms, 1);
isH = false(natoms, 1);

t(1) = tmp_t;

C = textscan(fid, '%s %s %d %f %f %f', natoms);
pos(C{3}, 1, 1) = C{4};
pos(C{3}, 2, 1) = C{5};
pos(C{3}, 3,1 ) = C{6};

for i=1:natoms
    residue_ind(C{3}(i)) = str2double(C{1}{i}(1:end-3));
end

for i=1:natoms
    if strcmp(C{2}{i}, 'CA')
        isCA(C{3}(i)) = true;
    end
end

for i=1:natoms
    if C{2}{i}(1) == 'H'
        isH(C{3}(i)) = true;
    end
end

C = textscan(fid, '%f', dim);
box(1, :) = C{1}';

%%

for i=2:nsteps
    fprintf('step = %d \n', i);
    
    C = textscan(fid,'%s %f',1,'delimiter','=');
    t(i) = C{2};

    C = textscan(fid, '%d', 1);

    C = textscan(fid, '%s %s %d %f %f %f', natoms);
    pos(C{3}, 1, i) = C{4};
    pos(C{3}, 2, i) = C{5};
    pos(C{3}, 3, i) = C{6};

    C = textscan(fid, '%f', dim);
    box(i, :) = C{1}';
end

%%
fclose(fid);

