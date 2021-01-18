import codecs
import os.path

import setuptools

setuptools.setup(
    name="maxwell_slice",
    version="0.1.0",
    author="Jeremy Magland and Ali Simons",
    author_email="",
    description="Visualizing outputs of 3D Maxwell solver",
    url="https://github.com/flatironinstitute/maxwell-slice",
    packages=setuptools.find_packages(), # this will automatically add a package for every subdirectory that has a __init__.py file (I think)
    include_package_data=True,
    scripts=[
    ],
    install_requires=[
        'numpy'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ]
)
