import pytest

import numpy as np

from pyzernike import ZernikeDescriptor

class TestCoeffcients:

    @pytest.mark.parametrize("dtype", (np.float32, np.float64))
    @pytest.mark.parametrize("order", (8, 20))
    def test_cuboid(self, order : int, dtype : type):
        arr = np.zeros((50,50,50), dtype = dtype)
        arr[15:25, 5:15, 35:45] = 1

        descriptor = ZernikeDescriptor.fit(data = arr, order = order)

        coefficients = descriptor.get_coefficients()

        base = np.load(f"tests/data/cuboid_order{order}.npy")

        diff = np.abs(np.subtract(coefficients, base)).sum()
        assert diff <= np.finfo(np.float32).resolution