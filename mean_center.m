function new_pos = mean_center(pos)

if ndims(pos) == 2
    [m, ~] = size(pos);
    new_pos = pos - repmat(mean(pos), m, 1);
elseif ndims(pos) == 3
    [m, n, p] = size(pos);
    new_pos = zeros(size(pos));
    for i=1:p
        new_pos(:,:,i) = pos(:,:,i) - repmat(mean(pos(:,:,i)), m, 1);
    end
else
    disp('shape of pos is invalid');
    return
end

