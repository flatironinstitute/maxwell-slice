import os
import shutil
from _temporarydirectory import TemporaryDirectory
from _shellscript import ShellScript
import numpy as np

def sphere_scat(nnx = 20, nny = 20, nnz = 20):
    # Creates a temporary) working directory (will get cleaned up even if exception is raised)
    # if remove == False, the temporary directory does not get removed
    with TemporaryDirectory(remove=False) as tmpdir:
        print(f'Using temporary dir:{tmpdir}')

        thisdir = os.path.dirname(os.path.realpath(__file__))

        # generate the input file
        # in future: we'll have a controls.dat.j2 and fill in the args
        shutil.copyfile(f'{thisdir}/../scratchspace/SphereScat/controls.dat.jinja2', f'{tmpdir}/controls.dat.jinja2')
        shutil.copyfile(f'{thisdir}/../scratchspace/SphereScat/spinsimp.dat', f'{tmpdir}/spinsimp.dat')

        controls = f'{tmpdir}/controls.dat.jinja2'
        # Read in the file
        with open(controls, 'r') as file :
            filedata = file.read()

            # Replace the target string
        filedata = filedata.replace('xdim', str(nnx))
        filedata = filedata.replace('ydim', str(nny))
        filedata = filedata.replace('zdim', str(nnz))

            # Write the file out again
        with open(controls, 'w') as file:
            file.write(filedata)
        print(controls)

        # Here is where we run the fortran program using a bash script
        s = ShellScript(f'''
        #!/bin/bash
        cd {tmpdir}

        {thisdir}/../scratchspace/SphereScat/int2 > output.txt
        ''')
        s.start()
        s.wait()

        print(f'nnx is {nnx}')
        print(f'nny is {nny}')
        print(f'nnz is {nnz}')

        data = np.fromfile(f'{tmpdir}/data.bin', dtype='>d')
        field = np.fromfile(f'{tmpdir}/field.bin', dtype='>d')
        data = data.reshape((5, nnx, nny, nnz))
        return data, field

if __name__ == '__main__':
    data, field = sphere_scat()
    print(data.shape)
    print(field.shape)