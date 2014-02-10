clear all
close all

%% load data
load rmsd_data
load all_atom_data

nsteps = length(t);

pos = reshape_pos(pos, nsteps);

%pos = pos(logical(1-isH'), :, :);
pos = pos(logical(isCA)', :, :);

pos = mean_center(pos);

rg = find_rg(pos);

%% PCA
pca_data = zeros(nsteps, size(pos,1) * size(pos,2));
template= pos(:, :, end);

for i=1:nsteps
    P = pos(:,:,i);
    Pnew = Kabsch(P, template);
    pca_data(i,:) = reshape(Pnew, 1, []);
end

pca_data = pca_data - repmat(mean(pca_data), nsteps, 1);

[V, D] = PCA(pca_data, 10);

for i=1:3
    figure; 
    plot(t, pca_data*V(:,i),'.')
end

for i=1:3
    figure; 
    plot(rg, pca_data*V(:,i),'.')
end

for i=1:3
    figure; 
    plot(rmsd, pca_data*V(:,i),'.')
end


return

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

[V, D] = dmaps(W, eps, 10);

%% plots
figure;
plot(diag(D),'.')

figure; 
plot(t, V(:,2),'.')

figure; 
plot(t, V(:,3),'.')

figure;
plot3(V(:,2),V(:,3),V(:,4),'.')


figure; 
plot(t, rg)

figure; 
plot(rg, V(:,2), '.')

figure; 
plot(rg, V(:,3), '.')

figure;
plot(rmsd, V(:,2),'.')

figure;
plot(rmsd, V(:,3),'.')
