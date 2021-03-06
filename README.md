# maxwell-slice
Placeholder for Ali/Jeremy/Leslie project for visualizing outputs of 3D Maxwell solver

## Possible tools

* K3D jupyter extension
    - [K3D jupyter extension](https://github.com/K3D-tools/K3D-jupyter) - tried it and it worked out of the box: seems pretty easy to 
    - Here's a project that uses k3d jupyter: https://discretisedfield.readthedocs.io/en/latest/ipynb/field-k3d-visualisation.html
    - Not sure how easy it will be to customize outside of the jupyter extension (custom react app)

* Vtk.js and ParaView
    - [react-vtkjs-viewport](https://github.com/OHIF/react-vtkjs-viewport) - I had some trouble compiling this into a react application
    - [paraview](https://www.paraview.org/)
    - [ParaViewWeb](https://www.paraview.org/web/)

* ipyvolume
    - [ipyvolume](https://github.com/maartenbreddels/ipyvolume)

* ipygany
    - [ipygany](https://github.com/QuantStack/ipygany)

* matplotlib
    - [mplot3d](https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html#d-plots-in-3d)

## Running SphereScat

Prerequisites: gfortran

```
cd scratchspace/SphereScat
./testcompile4
```

This will produce a bunch of .o files and an executable called `int2`

Let's run the program:

```
./int2
```

This reads input parameters from the following two files:
* controls.dat
* spinsimp.dat

And produces a few MATLAB files:

* data.m - 4 columns: x, y, z, scalar
    - n = number of rows = number of targets
    - represents phi
* field.m - 6 columns: re1, im1, re2, im2, re3, im3 (not reshaped)
    - n rows
    - represents grad(phi)
* sphereplot.m
    - sphere locations
* sphereplot2.m
    - sphere locations with phi plotted on the surface

![Output of sphereplot2](https://user-images.githubusercontent.com/3679296/104742094-06cf1b80-5718-11eb-9924-4e5ce70d929d.png)

*Output from sphereplot2.m*

## Install package in development (editable) mode

When installing in development mode, you don't need to reinstall every time you make changes to the source code.

```bash
# Clone this repo and cd to the base directory
pip install -e .
```

To test that it is installed. In ipython:

```
import maxwell_slice as ms
ms.test_function()
# should print "hello"
```

## How to test current sphere_scat program

At the moment, to run the sphere_scat program, all you need to do is:
```
python sphere_scat.py
```
Or (depending on how python is set up):
```
python3 sphere_scat.py
```
This should run the fortran code and print out a bunch of specifications from that code. It should also print out "test input arguments" at the very end.