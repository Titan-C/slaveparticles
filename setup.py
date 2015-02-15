# -*- coding: utf-8 -*-
"""
@author: Óscar Nájera
Installing packages on Slave Particles Mean Field Theories
"""
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import sys
import slaveparticles

class PyTest(TestCommand):
    """Test class to do test coverage analysis"""
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = ['--cov-report', 'term-missing',
                          '--cov', 'slaveparticles', 'tests/']
        self.test_suite = True
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


with open('README.rst') as f:
    long_description = f.read()


setup(
    name="slaveparticles",
    description="Educative code on Slave Particles",
    long_description=long_description,
    version=slaveparticles.__version__,
    packages=find_packages(),
    url="https://github.com/Titan-C/slaveparticles",
    author="Óscar Nájera",
    author_email='najera.oscar@gmail.com',
    license="GNU General Public License v3 (GPLv3)",

    install_requires=['numpy', 'scipy', 'matplotlib'],
    setup_requires=['Sphinx'],
    tests_require=['pytest', 'pytest-cov'],
    cmdclass={'test': PyTest},
)
