import ctypes
from typing import Sequence, Mapping, Dict, Tuple, Optional
from os import path
from sys import platform
from ctypes import cdll

import numpy as np
from numpy.ctypeslib import ndpointer

from mesh_to_depth import log

ctypes_lib_dirnames = [
    path.realpath(path.dirname(__file__)),
    path.realpath(path.join(path.dirname(__file__), '../../cpp/cmake-build-release/ctypes')),
]

lib = None


def find_lib_file(basename):
    for dirname in ctypes_lib_dirnames:
        candidate = path.join(dirname, basename)
        if path.isfile(candidate):
            return candidate
    log.error('Could not find {} in any of the following paths: \n{}'.format(basename, '\n'.join(ctypes_lib_dirnames)))


if platform.startswith('linux'):
    lib_filename = find_lib_file('libmesh2depth.so')
elif platform == "darwin":
    lib_filename = find_lib_file('libmesh2depth.dylib')
else:
    raise NotImplemented(platform)

try:
    lib = cdll.LoadLibrary(lib_filename)

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

    c_func = getattr(lib, 'mesh2depth_simple')
    c_func.argtypes = [
        ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),  # vertices. (num_vertices, 3)
        ctypes.c_size_t,  # num_vertices
        ndpointer(ctypes.c_uint32, flags="C_CONTIGUOUS"),  # faces. (num_faces, 3)
        ctypes.c_size_t,  # num_faces
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # cameras. (num_cameras, 16)
        ctypes.c_size_t,  # num_cameras
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # render params. (num_cameras, 3)
        ctypes.POINTER(ctypes.POINTER(ctypes.c_float)),  # output depth values. Dynamically allocated and needs to be freed.
    ]

except Exception as ex:
    import traceback

    traceback.print_exc()
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

    # copy() is important. Otherwise values will be lost after freeing memory.
    out_y_arr = [np.ctypeslib.as_array(out_y[i], shape=resolutions[i]).copy() for i in range(len(resolutions))]
    out_x_arr = [np.ctypeslib.as_array(out_x[i], shape=resolutions[i]).copy() for i in range(len(resolutions))]

    # Free dynamically allocated arrays.
    if cleanup:
        _free_arrays(out_y)
        _free_arrays(out_x)

    return out_y_arr, out_x_arr


def perspective_frustum(hw_ratio, x_fov, znear, zfar):
    """
    :param hw_ratio: height/width. e.g. 0.75
    :param x_fov: must be half-angle, in radians.
    :param znear: near clipping distance.
    :param zfar: far clipping distance.
    """
    assert znear != zfar
    right = np.abs(np.tan(x_fov) * znear)
    top = right * hw_ratio
    left = -right
    bottom = -top
    return [left, right, bottom, top, znear, zfar]


def mesh2depth_simple(vertices: np.ndarray, faces: np.ndarray, cam_and_render_params: Sequence[dict], empty_pixel_value=np.nan):
    """
    :param vertices: An array of shape (_, 3) and type float32.
    :param faces: An array of shape (_, 3) and type uint32.
    :param cam_and_render_params: A list of param dicts:
        cam_pos, cam_lookat, cam_up, x_fov, near, far, height, width, is_depth

        x_fov is the left-end to right-end angle in radians.
        is_depth indicates whether the output is a depth map or a ray displacement map. Depth map if true.

    :param empty_pixel_value: Background pixels where there is no depth will be filled with this value.
    """
    assert isinstance(vertices, np.ndarray)
    assert isinstance(faces, np.ndarray)
    assert isinstance(cam_and_render_params, (list, tuple))
    assert vertices.dtype == np.float32, 'dtype vertices must be np.float32. e.g. `vertices.astype(np.float32)`'
    assert faces.dtype == np.uint32, 'dtype faces must be np.uint32. e.g. `faces.astype(np.uint32)`'
    assert len(cam_and_render_params) > 0
    assert isinstance(cam_and_render_params[0], dict)
    assert vertices.flags['C_CONTIGUOUS'], 'vertices data must be contiguous in memory. e.g. See `numpy.ascontiguousarray` or do `vertices.copy()`'
    assert faces.flags['C_CONTIGUOUS'], 'faces data must be contiguous in memory. e.g. See `numpy.ascontiguousarray` or do `faces.copy()`'
    assert vertices.shape[1] == 3
    assert faces.shape[1] == 3

    c_func_name = 'mesh2depth_simple'
    num_cameras = len(cam_and_render_params)
    out_depth_values_ptr = (ctypes.POINTER(ctypes.c_float) * num_cameras)()

    cam_param_array_double = np.zeros((num_cameras, 16), dtype=np.float64, order='C')
    render_param_array_double = np.zeros((num_cameras, 3), dtype=np.float64, order='C')

    for i in range(num_cameras):
        cam_pos = cam_and_render_params[i]['cam_pos']
        cam_lookat = cam_and_render_params[i]['cam_lookat']
        cam_up = cam_and_render_params[i]['cam_up']
        x_fov = cam_and_render_params[i]['x_fov'] / 2  # The C++ implementation takes half angles.
        near = cam_and_render_params[i]['near']
        far = cam_and_render_params[i]['far']
        height = cam_and_render_params[i]['height']
        width = cam_and_render_params[i]['width']
        is_depth = bool(cam_and_render_params[i]['is_depth'])

        cam_param_array_double[i, 0:3] = cam_pos
        cam_param_array_double[i, 3:6] = cam_lookat
        cam_param_array_double[i, 6:9] = cam_up
        cam_param_array_double[i, 9:10] = 1.0  # is_perspective
        cam_param_array_double[i, 10:16] = perspective_frustum(hw_ratio=float(height) / width, x_fov=x_fov, znear=near, zfar=far)

        render_param_array_double[i, 0:2] = (height, width)
        render_param_array_double[i, 2] = float(is_depth)

    c_func = getattr(lib, c_func_name)

    c_func(
        vertices,
        int(vertices.shape[0]),
        faces,
        int(faces.shape[0]),
        cam_param_array_double,
        int(cam_param_array_double.shape[0]),
        render_param_array_double,
        out_depth_values_ptr,
    )

    out_depth_maps = []
    for i in range(num_cameras):
        height, width = int(render_param_array_double[i][0]), int(render_param_array_double[i][1])
        # copy() is important. Otherwise values will be lost after freeing memory.
        depth_map = np.ctypeslib.as_array(out_depth_values_ptr[i], shape=(height, width)).copy()
        if not np.isnan(empty_pixel_value):
            depth_map[np.isnan(depth_map)] = empty_pixel_value
        out_depth_maps.append(depth_map)

    _free_arrays(out_depth_values_ptr)

    return out_depth_maps
