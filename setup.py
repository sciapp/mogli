"""
Setup for mogli: https://pypi.python.org/pypi/mogli
"""
from setuptools import setup
import os.path
import shutil
import subprocess


def create_required_files():
    """
    Creates the files required for building a package.
    """
    # Manifest
    if not os.path.isfile('MANIFEST.in'):
        with open('MANIFEST.in', 'w') as manifest_in:
            manifest_in.write('include *.txt\n')
            manifest_in.write('include *.rst\n')

# License
if not os.path.isfile('LICENSE.txt'):
    shutil.copyfile('LICENSE', 'LICENSE.txt')

# Readme in reST format
if not os.path.isfile('README.rst'):
    subprocess.call(['pandoc', 'README.md', '-t', 'rst', '-o', 'README.rst'])

create_required_files()
# README.txt was created in the function above, so we can use its content for
# the long description:
with open('README.rst') as readme_file:
    long_description = readme_file.read()

setup(name='mogli',
      version='0.3.0',
      description='Simple visualization of molecules in python',
      long_description=long_description,
      author='Florian Rhiem',
      author_email='florian.rhiem@gmail.com',
      url='https://github.com/FlorianRhiem/mogli',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Visualization',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Physics',
          ],
      py_modules=['mogli'],
      install_requires=['glfw', 'gr', 'numpy', 'PyOpenGL'])
