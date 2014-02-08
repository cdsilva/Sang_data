function new_pos = reshape_pos(pos, nsteps)

[m, dim] = size(pos);
natoms = m / nsteps;

if mod(natoms, 1) ~= 0
    disp('nsteps is invalid')
    return
end

new_pos = reshape(pos, [natoms, dim, nsteps]);

