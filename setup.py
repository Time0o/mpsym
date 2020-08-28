#!/usr/bin/env python3

import os
import platform
import re
import subprocess
import sys

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
CMAKE_PROJECT_SETTINGS = os.path.join(SCRIPT_DIR, 'cmake', 'ProjectSettings.cmake')


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
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
            '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}'.format(extdir),
            '-DCMAKE_BUILD_TYPE=Release',
            '-DPYTHON_BINDINGS=ON',
            '-DPYTHON_EXECUTABLE=' + sys.executable,
            '-DEMBED_LUA=ON'
        ]

        subprocess.check_call(cmake_cmd, cwd=self.build_temp)

    def _build(self, ext):
        cmake_build_cmd = [
            'cmake',
            '--build', '.',
            '--config', 'Release',
            '--', '-j', str(os.cpu_count())
        ]

        subprocess.check_call(cmake_build_cmd, cwd=self.build_temp)


# check if CMake is available
try:
    subprocess.check_output(['cmake', '--version'])
except OSError:
    raise RuntimeError("cmake command must be available")


# parse project settings
project_settings = {}

r = re.compile(r'set\((.*) "(.*)"\)')

with open(CMAKE_PROJECT_SETTINGS, 'r') as f:
    for line in f:
        m = r.match(line.rstrip())

        var = m[1]
        var = var[len('CMAKE_PROJECT_'):]

        val = m[2]

        project_settings[var] = val


# setup
setup(
    name='py' + project_settings['NAME'],
    version=project_settings['VERSION'],
    description=project_settings['DESCRIPTION'],
    url=project_settings['HOMEPAGE_URL'],
    author=project_settings['AUTHOR'],
    license=project_settings['LICENSE'],
    cmdclass=dict(build_ext=CMakeBuild),
    ext_modules=[CMakeExtension(project_settings['NAME'])],
    zip_safe=False
)
