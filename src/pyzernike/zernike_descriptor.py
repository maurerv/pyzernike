#!python3
""" Implements 3D Zernike Descriptors.

    Copyright (c) 2024 European Molecular Biology Laboratory

    Author: Valentin Maurer <valentin.maurer@embl-hamburg.de>
"""

import numpy as np
from numpy.typing import NDArray


class ZernikeDescriptor:
    """
    Compute and interface with 3D Zernike Descriptors [1]_.

    References
    ----------
    ..[1] M. Novotni, R. Klein; Computer Aided Design 2004; 36(11):1047-1062

    """

    def __init__(self, descriptor: object):
        """
        Initialize the ZernikeDescriptor.

        Parameters
        ----------
        descriptor : object
            The underlying cpp descriptor object.
        """
        self._descriptor = descriptor

    @classmethod
    def fit(cls, data: NDArray, order: int = 20) -> "ZernikeDescriptor":
        """
        Compute the Zernike Descriptor for the given data.

        Parameters
        ----------
        data : NDArray
            The input 3D data array. Must be a cube with equal dimensions.
        order : int
            The order of the Zernike Descriptor. Must be a positive integer.
            20 is a value commonly found in the literature and thus set as default.

        Returns
        -------
        ZernikeDescriptor
            An instance of the ZernikeDescriptor class.

        Raises
        ------
        ValueError
            If the order is not a positive integer, if the input data is not a 3D cube,
            or if no data values larger than zero are found.
        NotImplementedError
            If the descriptor is not implemented for the given data type.
        """
        from ._core import ZernikeDescriptorFloat, ZernikeDescriptorDouble

        order = int(order)
        if order <= 0:
            raise ValueError("Order needs to be a positive integer.")

        if data.dtype == np.float32:
            func = ZernikeDescriptorFloat
        elif data.dtype == np.float64:
            func = ZernikeDescriptorDouble
        else:
            raise NotImplementedError(f"Descriptor not implemented for {data.dtype}.")

        if data.ndim != 3:
            raise ValueError(f"data needs to have ndim 3, got {data.ndim}.")

        if data.shape[1] != data.shape[0] or data.shape[2] != data.shape[0]:
            raise ValueError("Input must be a cube.")

        if np.sum(data >= 0) <= 0:
            raise ValueError("No data values larger than zero found.")

        descriptor = func(data, order=order)
        return cls(descriptor=descriptor)

    def get_coefficients(self) -> NDArray:
        """
        Return the computed Zernike descriptors

        Returns
        -------
        NDArray
            An array of Zernike Descriptor coefficients.
        """
        return self._descriptor.get_descriptors()

    def save_invariants(self, filepath: str) -> None:
        """
        Save the Zernike Descriptor invariants to a file.

        Parameters
        ----------
        filepath : str
            The path to the file where the invariants will be saved.
        """
        return self._descriptor.save_invariants(filepath)

    def reconstruct(self, box_size: int) -> NDArray:
        """
        Reconstruct 3D data from the Zernike Descriptor in a new box.

        Parameters
        ----------
        box_size : int
            The size of the reconstruction box.

        Returns
        -------
        NDArray
            The reconstructed 3D data array.
        """
        return self._descriptor.reconstruct(dim=box_size)
