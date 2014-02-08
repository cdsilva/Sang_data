function new_pos = mean_center(pos)

[m, n, p] = size(pos);

new_pos = pos;

for i=1:p
    new_pos(:,:,i) = new_pos(:,:,i) - repmat(mean(new_pos(:,:,i)), m, 1);
end

