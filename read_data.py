import numpy as np
import scipy.io

def read_coordinates(filename, out_mat_filename):
    dim = 3
    
    f = open(filename, 'r')

    t = []
    all_pos = []
    all_box = []
    
    line = f.readline()

    while line != '':
        t.append(float(line.split()[-1]))

        n = int(f.readline())
    
        pos = np.zeros((n, dim))
        atom_type = [[] for i in range(n)]
        
        for i in range(n):
            line = f.readline()
            ind = int(line.split()[2]) - 1
            atom_type[ind] = line.split()[1]
            for j in range(dim):
                pos[ind][j] = float(line.split()[3+j])

        all_pos.append(pos)

        line = f.readline()
        all_box.append([float(k) for k in line.split()])

        line = f.readline()
        
    f.close()

    isH = np.zeros(n)
    for i in range(n):
        if atom_type[i][0] == 'H':
            isH[i] = 1

    isCA = np.zeros(n)
    for i in range(n):
        if atom_type[i] == 'CA':
            isCA[i] = 1

    np.save('box', np.array(all_box))
    np.save('time', np.array(t))
    np.save('pos', np.array(all_pos))
    np.save('isH', isH)
    np.save('isCA', isCA)
    
    scipy.io.savemat(out_mat_filename,{'pos':np.array(all_pos[::10]), 't':np.array(t[::10]), 'box':np.array(all_box[::10]), 'isH':isH, 'isCA':isCA})

def read_rmsd(filename, output_filename):
    f = open(filename, 'r')

    for i in range(5):
        f.readline()

    t = []
    rmsd = []

    line = f.readline()
    while line != '':
        t.append(float(line.split()[0]))
        rmsd.append(float(line.split()[1]))
        line = f.readline()

    np.save('t_rmsd.npy', np.array(t))
    np.save('rmsd.npy', np.array(rmsd))

    scipy.io.savemat(output_filename,{'t':np.array(t[::10]), 'rmsd':np.array(rmsd[::10])})
    
if __name__ == '__main__':
    filename = 'all_atom_coordinate.gro'
    output_filename = 'all_atom_data.mat'

    rmsd_filename = 'rmsd_c_alpha.xvg'
    rmsd_output_filename = 'rmsd_data.mat'
    
    read_coordinates(filename, output_filename)
    read_rmsd(rmsd_filename, rmsd_output_filename)
