# 3D Zernike Descriptors

pyzernike is a Python library for computing 3D Zernike invariants as descriptors for 3D shape comparison.

## Usage

Here's a basic example of how to use the library:

```python
import numpy as np
from pyzernike import ZernikeDescriptor

# Create or load your 3D array here
arr = np.zeros((50,50,50), dtype = np.float32)
arr[15:25, 5:15, 35:45] = 1

# Fit the Zernike descriptor up to order 8
descriptor = ZernikeDescriptor.fit(data = arr, order = 8)

# Get Zernike coefficients
coefficients = descriptor.get_coefficients()

# Reconstruct the array
arr_rec = descriptor.reconstruct(box_size = 50)

# Write the coefficients to a binary file
descriptor.save_invariants("coefficients.inv")

```

## Installation

We recommend installation using one of the following methods

| Method   | Command                                                 |
|----------|---------------------------------------------------------|
| PyPi     | `pip install pyzernike`                                 |
| Source   | `pip install git+https://github.com/maurerv/pyzernike`  |


## Background

pyzernike provides Python bindings to C code written by Marcin Novotni, which was distributed under a GPL license and provided with the [paper](https://cg.cs.uni-bonn.de/backend/v1/files/publications/novotni-2004-shape.pdf):

M. Novotni, R. Klein "Shape Retrieval using 3D Zernike Descriptors" Computer Aided Design 2004; 36(11):1047-1062

A copy of that code serving as the basis for this project was obtained from https://github.com/codingforfun/ZernikeMoments. pyzernike includes a range of modifications to the original code base to improve performance but faithfully implements the original derivations.

