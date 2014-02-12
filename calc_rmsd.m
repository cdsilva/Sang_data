function rmsd = calc_rmsd(pos, template)

nsteps = size(pos, 3);

rmsd = zeros(nsteps, 1);

pos = mean_center(pos);
template = mean_center(template);

for i=1:nsteps
    Pnew = Kabsch(pos(:,:,i), template);
    rmsd(i) = sum((Pnew(:)-template(:)).^2);
end

