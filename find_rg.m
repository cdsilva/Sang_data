function rg = find_rg(pos)

[m, n, p] = size(pos);

pos = mean_center(pos);

rg = zeros(p, 1);
for i=1:p
    rg(i) = sum(sum((pos(:,:,i)).^2));
end