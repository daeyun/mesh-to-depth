import sys
from mesh_to_depth import ctypes_mesh2depth


def main():
    resolutions = [
        (1000, 1000),
        (500, 500),
    ]
    ctypes_mesh2depth._meshgrid_debug(resolutions, cleanup=True)


def main_no_op():
    # Used for generating valgrind suppression file.
    pass


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
        main()
        main()
        main()
    else:
        main_no_op()
