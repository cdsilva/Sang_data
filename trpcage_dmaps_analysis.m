clear all
close all

%% read & save data

% %read rmsd file that contains rmsd of alpha carbons
% [t, rmsd] = read_rmsd('rmsd_c_alpha.xvg');
% %save rmsd in mat file
% save('rmsd_data.mat', 't', 'rmsd');
% 
% %read in folded configuration
% [folded_pos, t, residue_ind, isH, isCA, box] = read_gro_file('original_folded.gro');
% %save folded configuration in mat file
% save('folded_structure.mat', 'folded_pos');
% 
% %read in trajectory
% [pos, t, residue_ind, isH, isCA, box] = read_gro_file('all_atom_coordinate.gro');
% %save trajectory in mat file
% save('all_atom_data.mat', 'pos','t', 'residue_ind', 'isH','isCA','box');

%% load saved data
% load rmsd_data
% 
% load folded_structure
% 
% load all_atom_data
%

%% subsample data & save

% %interval at which to subsample data
% step = 10;
% %start point for subsampling
% start = 1;
% 
% pos = pos(:,:,start:step:end);
% t = t(start:step:end);
% rmsd = rmsd(start:step:end);
% box = box(start:step:end,:);
% 
% save('all_atom_data_subsampled.mat', 'pos', 't', 'residue_ind', 'isH', 'isCA', 'box', 'folded_pos', 'rmsd');

%% load data
% 
load all_atom_data_subsampled

%% select which atoms to use

% calculate number of steps
nsteps = length(t);

% use all atoms except the hydrogens
ind = ~isH;
% use only alpha carbons
%ind = isCA;

% extract data pertaining to relevant atoms
pos = pos(ind, :, :);
folded_pos = folded_pos(ind, :, :);
residue_ind = residue_ind(ind);

% mean-center configurations
pos = mean_center(pos);

%% calculate different (potential) order parameters

% calculate rmsd of alpha-helix
alpha_helix_ind = find(residue_ind <= 8 & residue_ind >= 2);
rmsd_alpha_helix = calc_rmsd(pos(alpha_helix_ind,:,:), folded_pos(alpha_helix_ind, :));

% calculate rmsd of 3-10 helix
helix310_ind = find(residue_ind <= 14 & residue_ind >= 11);
rmsd_helix310 = calc_rmsd(pos(helix310_ind,:,:), folded_pos(helix310_ind, :));

% calculate length of salt bridge
salt_bridge = zeros(nsteps, 1);
ind16 = find(residue_ind == 16);
ind9 = find(residue_ind == 9);
for i=1:nsteps
    salt_bridge(i) = norm(mean(pos(ind16,:,i)) - mean(pos(ind9,:,i)));
end

%% construct distance matrix for DMAPS
W = zeros(nsteps);
h = waitbar(0, 'Computing pairwise distance matrix....');
for i=1:nsteps
    waitbar(i^2/nsteps^2, h);
    for j=1:i-1
        P = pos(:,:,i);
        Q = pos(:,:,j);
        
        % align configurations
        Pnew = Kabsch(P, Q);
        
        % calculate distance between aligned configurations
        W(i,j) = sum((Pnew(:)-Q(:)).^2);
        W(j,i) = W(i,j);
    end
end
close(h);

%% compute dmaps embedding

% use median of distances as epsilon for distance kernel
eps = median(W(:));
% compute dmaps embedding coordinates
% the coordinates are stored in the matrix dmaps_coords
% each row of the matrix corresponds to a data point, each column
% corresponds to a specific coordinate
% the first coordinate (the first column of dmaps_coords) is trivial (it
% will always be a constant vector), and so we don't use it
[dmaps_coords, D_dmaps] = dmaps(W, eps, 10);

%% plots

% plot eigenvalues
figure;
plot(diag(D_dmaps), '.')
xlabel('k')
ylabel('\lambda_k')
title('DMAPS')

% plot 2-D dmaps embedding, color by rmsd
figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,3),200,rmsd,'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('DMAPS: colored by rmsd')

% plot 2-D dmaps embedding, color by third (nontrivial) dmaps coordinate
figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,3),200,dmaps_coords(:,4),'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('DMAPS: colored by \phi_4')

% plot 3-D dmaps embedding, color by rmsd
figure;
scatter3(dmaps_coords(:,2), dmaps_coords(:,3), dmaps_coords(:,4),200,rmsd,'.')
xlabel('\phi_2')
ylabel('\phi_3')
zlabel('\phi_4')
title('DMAPS: colored by rmsd')

% plot 2-D dmaps embedding (in 1st and 3rd coordinate), color by rmsd
figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,4),200,rmsd,'.')
xlabel('\phi_2')
ylabel('\phi_4')
title('DMAPS: colored by rmsd')

% plot 2-D dmaps embedding, color by rmsd of alpha helix
figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,3),200,rmsd_alpha_helix,'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('DMAPS: colored by rmsd of alpha helix')

% plot 2-D dmaps embedding (in 1st and 3rd coordinate), color by rmsd of
% alpha helix
figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,4),200,rmsd_alpha_helix,'.')
xlabel('\phi_2')
ylabel('\phi_4')
title('DMAPS: colored by rmsd of alpha helix')

%% make figure with some points highlighted
figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,3),200,rmsd,'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('DMAPS: colored by rmsd')
hold on
ind = [100, 400, 1580, 2000, 2557, 3000];
plot(dmaps_coords(ind,2),dmaps_coords(ind,3),'.k','markersize',40)
for i=1:length(ind)
    text(dmaps_coords(ind(i),2)+0.002,dmaps_coords(ind(i),3)+0.002, num2str(t(ind(i))))
end
