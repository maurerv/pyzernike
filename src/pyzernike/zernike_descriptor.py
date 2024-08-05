import numpy as np

from numpy.typing import NDArray

class ZernikeDescriptor:
    def __init__(self, descriptor : object):
        self._descriptor = descriptor

    def get_coefficients(self):
        return self._descriptor.get_descriptors()

    def reconstruct(self, box_size : int) -> NDArray:
        return self._descriptor.reconstruct(dim = box_size)

    @classmethod
    def fit(cls, data : NDArray, order : int) -> "ZernikeDescriptor":
        from ._core import ZernikeDescriptorFloat, ZernikeDescriptorDouble

        if data.dtype == np.float32:
            func = ZernikeDescriptorFloat
        elif data.dtype == np.float64:
            func = ZernikeDescriptorDouble
        else:
            raise NotImplementedError(
                f"Descriptor not implemented for {data.dtype}."
            )

        descriptor = func(data, order = order)
        return cls(descriptor = descriptor)

