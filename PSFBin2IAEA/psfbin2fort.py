import numpy as np
from scipy.io import FortranFile
import argparse


#function to read binary and save as fortran routine
def fortcreate(inp, out, x0, y0, z0, ds):

    #read and reshape data
    vec = np.fromfile(inp, dtype=np.float32).reshape((-1,7))

    #place voxel center at origin
    vec[:,0] = vec[:,0] - x0 
    vec[:,1] = vec[:,1] - x0 
    vec[:,2] = vec[:,2] - x0

    #travel backwards
    vec[:,0] = vec[:,0]  - ds*vec[:,3]
    vec[:,1] = vec[:,1]  - ds*vec[:,4]
    vec[:,2] = vec[:,2]  - ds*vec[:,5]

    f = FortranFile(out, 'w')
    f.write_record(vec.T)
    f.close()

    nparticles = vec.shape[0]

    return nparticles
    
def arg_parser():
    parser = argparse.ArgumentParser(description='Converts PSF binary file to Fortran format')
    parser.add_argument('--i', type=str, nargs='+',
                        help='path to binary file', dest='inp', required=True)
    parser.add_argument('--o', metavar='O', type=str, nargs='+',
                        help='output Fortran file', dest='out', required=True)
    parser.add_argument('--x', metavar='X', type=float, nargs='+',
                        help='X voxel center position (in cm)', dest='x0', default=[0.0], required=False)
    parser.add_argument('--y', metavar='Y', type=float, nargs='+',
                        help='Y voxel center position (in cm)', dest='y0', default=[0.0], required=False)
    parser.add_argument('--z', metavar='Z', type=float, nargs='+',
                        help='Z voxel center position (in cm)', dest='z0', default=[0.0], required=False)
    parser.add_argument('--s', metavar='S', type=float, nargs='+',
                        help='Backwards step (in cm)', dest='ds', default=[0.0], required=False)
    return parser


if __name__ == "__main__":
    parser = arg_parser().parse_args()


    nparticles = fortcreate(parser.inp[0], parser.out[0], parser.x0[0], parser.y0[0], parser.z0[0], parser.ds[0])

    print('PSF Fortran file generated sucessfully with {} photons'.format(nparticles))
    print('File saved as {}'.format(parser.out[0]))

    
    
    




