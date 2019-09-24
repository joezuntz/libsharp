from setuptools import setup, Extension, Distribution
import setuptools.command.build_ext

import sys
import sysconfig
import distutils.sysconfig


# FIXME this has to be set from outside!
sharp_libpath='/home/martin/codes/sharp2/.libs'

class _deferred_pybind11_include(object):
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


def _remove_strict_prototype_option_from_distutils_config():
    strict_prototypes = '-Wstrict-prototypes'
    config = distutils.sysconfig.get_config_vars()
    for key, value in config.items():
        if strict_prototypes in str(value):
            config[key] = config[key].replace(strict_prototypes, '')


_remove_strict_prototype_option_from_distutils_config()


extra_cc_compile_args = []
include_dirs = ['../',
                _deferred_pybind11_include(),
                _deferred_pybind11_include(True)]

python_module_link_args = []

if sys.platform == 'darwin':
    extra_cc_compile_args += ['--std=c++11', '--stdlib=libc++',
                              '-mmacosx-version-min=10.9']
    vars = distutils.sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '')
    python_module_link_args += ['-bundle']
else:
    extra_cc_compile_args += ['--std=c++11']
    python_module_link_args += ["-Wl,-rpath,$ORIGIN"]


def get_extension_modules():
    return [Extension('pysharp',
                      sources=['pysharp.cc'],
                      include_dirs=include_dirs,
                      extra_compile_args=extra_cc_compile_args,
                      libraries=["sharp"],
                      library_dirs=[sharp_libpath],
                      extra_link_args=python_module_link_args)]


setup(name='pysharp',
      version='0.0.1',
      description='Python bindings for libsharp',
      include_package_data=True,
      author='Martin Reinecke',
      author_email='martin@mpa-garching.mpg.de',
      packages=[],
      setup_requires=['numpy>=1.10.4', 'pybind11>=2.2.1'],
      ext_modules=get_extension_modules(),
      install_requires=['numpy>=1.10.4', 'pybind11>=2.2.1']
      )
