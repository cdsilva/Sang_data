import numpy as np
import scipy.io

mat_stride = 4

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
        residue_ind = [[] for i in range(n)]
        
        for i in range(n):
            line = f.readline()
            ind = int(line.split()[2]) - 1
            atom_type[ind] = line.split()[1]
            residue_ind[ind] = int(line.split()[0][0:-3])
            
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
    np.save('residue_ind', np.array(residue_ind))
            
    scipy.io.savemat(out_mat_filename,{'pos':np.array(all_pos[::mat_stride]), 't':np.array(t[::mat_stride]), 'box':np.array(all_box[::mat_stride]), 'isH':isH, 'isCA':isCA, 'residue_ind':residue_ind})

def read_folded():

    filename = 'original_folded.gro'
    
    dim = 3
    
    f = open(filename, 'r')

    line = f.readline()
    n = int(f.readline())
    
    pos = np.zeros((n, dim))
    for i in range(n):
        line = f.readline()
        ind = int(line.split()[2]) - 1
                    
        for j in range(dim):
            pos[ind][j] = float(line.split()[3+j])

    f.close()

    np.save('folded_structure', pos)
    scipy.io.savemat('folded_structure.mat',{'folded_pos':pos})

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

    scipy.io.savemat(output_filename,{'t':np.array(t[::mat_stride]), 'rmsd':np.array(rmsd[::mat_stride])})
    
if __name__ == '__main__':
    filename = 'all_atom_coordinate.gro'
    output_filename = 'all_atom_data.mat'

    rmsd_filename = 'rmsd_c_alpha.xvg'
    rmsd_output_filename = 'rmsd_data.mat'
    
    read_coordinates(filename, output_filename)
    read_folded()
    read_rmsd(rmsd_filename, rmsd_output_filename)
