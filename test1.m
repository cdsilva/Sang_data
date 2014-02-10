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
load rmsd_data

load folded_structure
folded_pos = pos;

load all_atom_data

%%
step = 10;
start = 5000;
pos = pos(:,:,start:step:end);
t = t(start:step:end);
rmsd = rmsd(start:step:end);
box = box(start:step:end,:);

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

helix310_ind = find(residue_ind <= 14 & residue_ind >= 11);
rg_helix310 = find_rg(pos(helix310_ind,:,:));

salt_bridge = zeros(nsteps, 1);
ind16 = find(residue_ind == 16);
ind9 = find(residue_ind == 9);
for i=1:nsteps
    salt_bridge(i) = norm(mean(pos(ind16,:,i)) - mean(pos(ind9,:,i)));
end
%% calculate rmsd

natoms = size(pos, 1);

template = folded_pos - repmat(mean(folded_pos), natoms, 1);
rmsd_tmp = zeros(nsteps, 1);
for i=1:nsteps
    Pnew = Kabsch(pos(:,:,i), template);
    rmsd_tmp(i) = sum((Pnew(:)-template(:)).^2)/natoms;
end

figure;
plot(t, rmsd_tmp)

%% PCA
pca_data = zeros(nsteps, size(pos,1) * size(pos,2));

for i=1:nsteps
    P = pos(:,:,i);
    Pnew = Kabsch(P, template);
    pca_data(i,:) = reshape(Pnew, 1, []);
end
mean_data = mean(pca_data);

pca_data = pca_data - repmat(mean_data, nsteps, 1);

[V, D] = PCA(pca_data, 10);

figure;
plot(diag(D), '.')
xlabel('k')
ylabel('\lambda_k')
title('PCA')

for i=1:3
    figure;
    plot(t, pca_data*V(:,i),'.')
    xlabel('t')
    ylabel(sprintf('Projection onto PC %d', i))
    title('PCA')
end

figure;
scatter(pca_data*V(:,1), pca_data*V(:,2),200,t,'.')
xlabel('Projection onto PC 1')
ylabel('Projection onto PC 2')
title('PCA')

figure;
scatter3(pca_data*V(:,1), pca_data*V(:,2),  pca_data*V(:,3),200,t,'.')
xlabel('Projection onto PC 1')
ylabel('Projection onto PC 2')
zlabel('Projection onto PC 3')
title('PCA')

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
[V, D] = dmaps(W, eps, 10, 1e-5);

%% plots
figure;
plot(diag(D),'.')

for i=2:3
    figure;
    plot(t, V(:,i),'.')
    xlabel('t')
    ylabel(sprintf('\\phi_%d', i))
    title('DMAPS')
end

figure;
scatter3(V(:,2),V(:,3),V(:,4),200, t, '.')
xlabel('\phi_2')
ylabel('\phi_3')
zlabel('\phi_4')

for i=2:3
    figure;
    plot(rg, V(:,i), '.')
    xlabel('R_g')
    ylabel(sprintf('\\phi_%d', i))
    title('DMAPS')
end

for i=2:3
    figure;
    plot(rmsd, V(:,i),'.')
    xlabel('rmsd (folded structure)')
    ylabel(sprintf('\\phi_%d', i))
    title('DMAPS')
end

