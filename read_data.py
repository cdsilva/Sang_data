import numpy as np
import scipy.io

def read_file(filename, out_mat_filename):
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
    
    np.save('box', np.array(all_box))
    np.save('time', np.array(t))
    np.save('pos', np.array(all_pos))
    np.save('isH', isH)
    
    scipy.io.savemat(out_mat_filename,{'pos':np.array(all_pos[::10]), 't':np.array(t[::10]), 'box':np.array(all_box[::10]), 'isH':isH})

if __name__ == '__main__':
    filename = 'all_atom_coordinate.gro'
    output_filename = 'all_atom_data.mat'
    
    read_file(filename, output_filename)

