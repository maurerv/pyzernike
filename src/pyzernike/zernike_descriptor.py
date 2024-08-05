import numpy as np

from numpy.typing import NDArray


class ZernikeDescriptor:
    def __init__(self, descriptor: object):
        self._descriptor = descriptor

    @classmethod
    def fit(cls, data: NDArray, order: int) -> "ZernikeDescriptor":
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

        descriptor = func(data, order=order)
        return cls(descriptor=descriptor)

    def get_coefficients(self):
        return self._descriptor.get_descriptors()

    def save_invariants(self, filepath: str) -> None:
        return self._descriptor.save_invariants(filepath)

    def reconstruct(self, box_size: int) -> NDArray:
        return self._descriptor.reconstruct(dim=box_size)
