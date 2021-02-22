import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext


# Based on https://github.com/pybind/cmake_example/blob/master/setup.py
#          https://www.benjack.io/2018/02/02/python-cpp-revisited.html
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                               out.decode()).group(1))
        if cmake_version < '3.1.3':
            raise RuntimeError("CMake >= 3.1.3 is required")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        os.system('ls -lah /tmp/pip-req-build-*/')

        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if not self.debug:
            cmake_args.extend([
                '-Dtest=OFF',
                '-DDEBUG=OFF',
            ])

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
                cfg.upper(),
                extdir)]
            if sys.maxsize > 2 ** 32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp)
        print()  # Add an empty line for cleaner output


_email = lambda x: x.replace(' [at] ', '@')
_proj_root = os.path.dirname(os.path.realpath(__file__))

test_deps = [
    'pytest',
]

extras = {
    'test': test_deps,
}

setup(
    name='mesh-to-depth',
    version='0.1.2',
    author='Daeyun Shin',
    author_email=_email('daeyuns [at] uci.edu'),
    description='Generate depth maps, given a mesh and camera parameters',
    long_description='',
    url='https://github.com/daeyun/mesh_to_depth',
    packages=['mesh_to_depth'],
    package_dir={'': 'python'},
    ext_modules=[CMakeExtension('mesh_to_depth/cpp_lib', sourcedir=os.path.join(_proj_root, 'cpp'))],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'setuptools',
        'networkx',
        'trimesh',
        'pillow',
        # 'pyassimp>=4.0',
    ],
    tests_require=test_deps,
    extras_require=extras,
)
