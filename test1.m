clear all
close all

%% read data

% [t, rmsd] = read_rmsd('rmsd_c_alpha.xvg');
% save('rmsd_data.mat', 't', 'rmsd');
% 
% [pos, t, residue_ind, isH, isCA, box] = read_gro_file('original_folded.gro');
% save('folded_structure.mat', 'pos');
% 
% [pos, t, residue_ind, isH, isCA, box] = read_gro_file('all_atom_coordinate.gro');
% save('all_atom_data.mat', 'pos','t', 'residue_ind', 'isH','isCA','box');

%% load data
% load rmsd_data
% 
% load folded_structure
% folded_pos = pos;
% 
% load all_atom_data
% 
% step = 10;
% start = 1;
% pos = pos(:,:,start:step:end);
% t = t(start:step:end);
% rmsd = rmsd(start:step:end);
% box = box(start:step:end,:);
% 
% save('all_atom_data_subsampled.mat', 'pos', 't', 'residue_ind', 'isH', 'isCA', 'box', 'folded_pos', 'rmsd');

%%
load all_atom_data_subsampled

% step = 1;
% start = 500;
% pos = pos(:,:,start:step:end);
% t = t(start:step:end);
% rmsd = rmsd(start:step:end);
% box = box(start:step:end,:);

nsteps = length(t);

ind = ~isH;
%ind = isCA;

pos = pos(ind, :, :);
folded_pos = folded_pos(ind, :, :);
residue_ind = residue_ind(ind);

pos = mean_center(pos);

rg = find_rg(pos);

alpha_helix_ind = find(residue_ind <= 8 & residue_ind >= 2);
rg_alpha_helix = find_rg(pos(alpha_helix_ind,:,:));
rmsd_alpha_helix = calc_rmsd(pos(alpha_helix_ind,:,:), folded_pos(alpha_helix_ind, :));

helix310_ind = find(residue_ind <= 14 & residue_ind >= 11);
rg_helix310 = find_rg(pos(helix310_ind,:,:));
rmsd_helix310 = calc_rmsd(pos(helix310_ind,:,:), folded_pos(helix310_ind, :));

salt_bridge = zeros(nsteps, 1);
ind16 = find(residue_ind == 16);
ind9 = find(residue_ind == 9);
for i=1:nsteps
    salt_bridge(i) = norm(mean(pos(ind16,:,i)) - mean(pos(ind9,:,i)));
end

%% PCA
pca_data = zeros(nsteps, size(pos,1) * size(pos,2));

template = mean_center(folded_pos);
for i=1:nsteps
    P = pos(:,:,i);
    Pnew = Kabsch(P, template);
    pca_data(i,:) = reshape(Pnew, 1, []);
end
mean_data = mean(pca_data);

pca_data = mean_center(pca_data);

[V, D_pca] = PCA(pca_data, 10);
pca_coords = pca_data * V;

%% PCA plots
figure;
plot(diag(D_pca), '.')
xlabel('k')
ylabel('\lambda_k')
title('PCA')

figure;
scatter(pca_coords(:,1),pca_coords(:,2),200,rmsd,'.')
xlabel('Projection onto PC 1')
ylabel('Projection onto PC 2')
title('PCA: colored by rmsd')

figure;
scatter3(pca_coords(:,1), pca_coords(:,2), pca_coords(:,3),200,rmsd,'.')
xlabel('Projection onto PC 1')
ylabel('Projection onto PC 2')
zlabel('Projection onto PC 3')
title('PCA: colored by rmsd')

figure;
scatter3(pca_coords(:,1), pca_coords(:,2), pca_coords(:,3),200,rmsd_helix310,'.')
xlabel('Projection onto PC 1')
ylabel('Projection onto PC 2')
zlabel('Projection onto PC 3')
title('PCA: colored by rmsd of 3-10 helix')


%% play movie

% figure;
% for i=1:nsteps
%     plot3(pos(:,1,i), pos(:,2,i), pos(:,3,i))
%     axis([-2 2 -2 2 -2 2])
%     pause(0.1)
% end

%% construct distance matrix
W = zeros(nsteps);
for i=1:nsteps
    i
    for j=1:i-1
        P = pos(:,:,i);
        Q = pos(:,:,j);
        Pnew = Kabsch(P, Q);
        
        W(i,j) = sum((Pnew(:)-Q(:)).^2);
        W(j,i) = W(i,j);
    end
end

%% dmaps

eps = median(W(:));
[dmaps_coords, D_dmaps] = dmaps(W, eps, 10);

%% plots
figure;
plot(diag(D_dmaps), '.')
xlabel('k')
ylabel('\lambda_k')
title('DMAPS')



figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,3),200,rmsd,'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('DMAPS: colored by rmsd')

figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,3),200,dmaps_coords(:,4),'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('DMAPS: colored by \phi_4')

figure;
scatter3(dmaps_coords(:,2), dmaps_coords(:,3), dmaps_coords(:,4),200,rmsd,'.')
xlabel('\phi_2')
ylabel('\phi_3')
zlabel('\phi_4')
title('DMAPS: colored by rmsd')

figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,4),200,rmsd,'.')
xlabel('\phi_2')
ylabel('\phi_4')
title('DMAPS: colored by rmsd')

figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,3),200,rmsd_alpha_helix,'.')
xlabel('\phi_2')
ylabel('\phi_3')
title('DMAPS: colored by rmsd of alpha helix')

figure;
scatter(dmaps_coords(:,2),dmaps_coords(:,4),200,rmsd_alpha_helix,'.')
xlabel('\phi_2')
ylabel('\phi_4')
title('DMAPS: colored by rmsd of alpha helix')
