#!/usr/bin/env python3

import os
import subprocess
import sys

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    user_options = build_ext.user_options + [
        ('debug-build', None, "Build in debug mode"),
        ('cmake-extra-opts=', None, "Additional options passed to CMake"),
    ]

    def initialize_options(self):
        build_ext.initialize_options(self)
        self.debug_build = False
        self.cmake_extra_opts = None

    def finalize_options(self):
        build_ext.finalize_options(self)
        self.debug_build = bool(self.debug_build)

    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("cmake command must be available")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        self._setup(ext)
        self._build(ext)

    def _setup(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_cmd = [
            'cmake',
            ext.sourcedir,
            '-DCMAKE_BUILD_TYPE=' + ('RelWithDebInfo' if self.debug_build else 'Release'),
            '-DIGNORE_INSTALL_PREFIX=ON',
            '-DLINK_STATIC=ON',
            '-DPYTHON_BINDINGS=ON',
            '-DPYTHON_MODULE_INTERNAL_OUT_DIR={}'.format(extdir),
            '-DPYTHON_EXECUTABLE=' + sys.executable,
            '-DLUA_EMBED=ON'
        ]

        if self.cmake_extra_opts is not None:
            cmake_cmd += self.cmake_extra_opts.split()

        subprocess.check_call(cmake_cmd, cwd=self.build_temp)

    def _build(self, ext):
        cmake_build_cmd = [
            'cmake',
            '--build', '.',
            '--config', 'Release',
            '--', '-j', str(os.cpu_count())
        ]

        subprocess.check_call(cmake_build_cmd, cwd=self.build_temp)


# Hack that enables us to run both 'python setup.py install' and 'pip install .'
# from this directory. Without this, we would not be able to locate this repo's
# root directory from the temporary build directory created by pip.

builddir = os.path.basename(os.path.dirname(sys.argv[0]))

if builddir.startswith('pip-req-build'):
    if os.environ.get('TMPDIR') != '.':
        raise ValueError("For technical reasonds, TMPDIR must be set to '.' when building with pip")

    sourcedir = '../..'
else:
    sourcedir = '..'


setup(
    name='mpsym',
    version='0.5',
    description="MPSoC Symmetry Reduction",
    long_description="mpsym is a C++/Lua/Python library that makes it possible to determine whether mappings of computational tasks to multiprocessor systems are equivalent by symmetry. It can also potentially be used to solve more general graph symmetry problems.",
    url="https://github.com/Time0o/mpsym",
    author="Timo Nicolai",
    author_email="timonicolai@arcor.de",
    license="MIT",
    cmdclass=dict(build_ext=CMakeBuild),
    ext_modules=[CMakeExtension('mpsym._mpsym', sourcedir)],
    packages=['mpsym'],
    setup_requires=['wheel'],
    zip_safe=False
)
