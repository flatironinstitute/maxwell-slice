import os
import shutil
from ._temporarydirectory import TemporaryDirectory
from ._shellscript import ShellScript
from jinja2 import Template
import numpy as np

def sphere_scat(nnx = 21, nny = 22, nnz = 23):
    # Creates a temporary) working directory (will get cleaned up even if exception is raised)
    # if remove == False, the temporary directory does not get removed
    with TemporaryDirectory(remove=False) as tmpdir:
        print(f'Using temporary dir:{tmpdir}')

        thisdir = os.path.dirname(os.path.realpath(__file__))

        # generate the input file
        # in future: we'll have a controls.dat.j2 and fill in the args
        shutil.copyfile(f'{thisdir}/../scratchspace/SphereScat/spinsimp.dat', f'{tmpdir}/spinsimp.dat')

        with open(f'{thisdir}/../scratchspace/SphereScat/controls.dat.jinja2', 'r') as f:
            controls_template = Template(f.read())
        with open(f'{tmpdir}/controls.dat', 'w') as f:
            f.write(controls_template.render(nnx=nnx, nny=nny, nnz=nnz))
        
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
        data = data.reshape((7, nnx, nny, nnz), order='F')
        field = field.reshape((6, nnx, nny, nnz), order='F')
        field = field[[0, 2, 4]] + 1j * field[[1, 3, 5]]
        return data, field