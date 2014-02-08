clear all
close all

load all_atom_data

nsteps = length(time);

pos = reshape_pos(pos, nsteps);

pos = mean_center(pos);

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

