import numpy as np
from mesh_to_depth import ctypes_mesh2depth


def test_meshgrid():
    resolutions = [
        (12, 16),
        (24, 36),
        (1000, 1000),
        (2, 3),
        (1, 1),
    ]

    x1_list, x2_list = ctypes_mesh2depth._meshgrid_debug(resolutions, cleanup=True)

    for i in range(len(resolutions)):
        x1, x2 = np.meshgrid(np.arange(resolutions[i][0]), np.arange(resolutions[i][1]), indexing='ij')
        assert x1_list[i].shape == resolutions[i], i
        assert x2_list[i].shape == resolutions[i], i
        assert np.all(x1_list[i] == x1), i
        assert np.all(x2_list[i] == x2), i
