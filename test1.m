clear all
close all

%%
load all_atom_data

nsteps = length(t);

pos = reshape_pos(pos, nsteps);

pos = pos(logical(1-isH'), :, :);
%pos = pos(logical(isCA)', :, :);

pos = mean_center(pos);

%%
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

%%
eps = 0.25 * median(W(:));

[V, D] = dmaps(W, eps, 10);

%%
figure;
plot(diag(D),'.')

figure; 
plot(t, V(:,2),'.')

figure; 
plot(t, V(:,3),'.')

figure;
plot3(V(:,2),V(:,3),V(:,4),'.')

rg = find_rg(pos);

figure; 
plot(t, rg)

figure; 
plot(rg, V(:,2), '.')

figure; 
plot(rg, V(:,3), '.')

