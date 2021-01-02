import ctypes
from typing import Sequence, Mapping, Dict, Tuple, Optional
from os import path
from sys import platform
from ctypes import cdll

import numpy as np
from numpy.ctypeslib import ndpointer

from mesh_to_depth import log

ctypes_lib_dirname = path.realpath(path.join(path.dirname(__file__), '../../cpp/cmake-build-release/ctypes'))

lib = None
if platform.startswith('linux'):
    lib_filename = path.join(ctypes_lib_dirname, 'libmesh2depth.so')
elif platform == "darwin":
    lib_filename = path.join(ctypes_lib_dirname, 'libmesh2depth.dylib')
else:
    raise NotImplemented(platform)

if path.isfile(lib_filename):
    lib = cdll.LoadLibrary(lib_filename)
else:
    log.error('file does not exist: {}'.format(lib_filename))
if lib:
    c_func = getattr(lib, '_meshgrid_debug')
    c_func.argtypes = [
        ctypes.c_size_t,  # num_images
        ndpointer(ctypes.c_uint32, flags="C_CONTIGUOUS"),  # resolutions. (_, 2)
        ctypes.POINTER(ctypes.POINTER(ctypes.c_uint32)),  # out_y
        ctypes.POINTER(ctypes.POINTER(ctypes.c_uint32)),  # out_x
    ]

    c_func = getattr(lib, 'free_arrays')
    c_func.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),  # arr_ptr. Array of pointers.
        ctypes.c_size_t,  # size
    ]
else:
    log.error('External library not loaded correctly: {}'.format(lib_filename))


def is_available():
    return lib is not None


def _free_arrays(arr_ptr: ctypes.Array):
    isinstance(arr_ptr, ctypes.Array)

    num_arrays = len(arr_ptr)
    arr_ptr_void = ctypes.cast(arr_ptr, ctypes.POINTER(ctypes.c_void_p))

    c_func_name = 'free_arrays'
    c_func = getattr(lib, c_func_name)
    c_func(
        arr_ptr_void,
        num_arrays,
    )


def _meshgrid_debug(resolutions: list, cleanup=True) -> Tuple[list, list]:
    """
    :param resolutions: List of (h, w) image sizes.
    :param cleanup: Set to False to cause memory leak.
    :return: (y, x). Lists of grid coordinates.
    """
    assert isinstance(resolutions, (list, tuple))
    assert len(resolutions) > 0
    assert len(resolutions[0]) == 2

    c_func_name = '_meshgrid_debug'
    num_images = len(resolutions)
    resolutions_arr = np.array(resolutions, dtype=np.uint32)
    out_y = (ctypes.POINTER(ctypes.c_uint32) * num_images)()
    out_x = (ctypes.POINTER(ctypes.c_uint32) * num_images)()

    c_func = getattr(lib, c_func_name)

    c_func(
        num_images,
        resolutions_arr,
        out_y,
        out_x,
    )

    # copy() is important. Otherwise values will be lose after freeing memory.
    out_y_arr = [np.ctypeslib.as_array(out_y[i], shape=resolutions[i]).copy() for i in range(len(resolutions))]
    out_x_arr = [np.ctypeslib.as_array(out_x[i], shape=resolutions[i]).copy() for i in range(len(resolutions))]

    # Free dynamically allocated arrays.
    if cleanup:
        _free_arrays(out_y)
        _free_arrays(out_x)

    return out_y_arr, out_x_arr
