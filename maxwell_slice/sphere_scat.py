from ._temporarydirectory import TemporaryDirectory
from ._shellscript import ShellScript

def sphere_scat():
    # Creates a temporary) working directory (will get cleaned up even if exception is raised)
    # if remove == False, the temporary directory does not get removed
    with TemporaryDirectory(remove=True) as tmpdir:
        input_fname = tmpdir + '/input.txt'
        output_fname = tmpdir + '/output.txt'

        # generate input files in working directory
        with open(input_fname, 'w') as f:
            f.write('test input arguments')
        # run the fortran program

        # Here is where we run the fortran program using a bash script
        s = ShellScript(f'''
        #!/bin/bash

        cd {tmpdir}
        # this would be replaced by the fortran call...
        cat input.txt > output.txt
        ''')
        s.start()
        s.wait()

        # read output
        with open(output_fname, 'r') as f:
            x = f.read()
        # return the output
        return x