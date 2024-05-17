import subprocess

from setuptools import find_packages, setup
from setuptools.command.install import install


class CustomInstallCommand(install):
    def run(self):
        subprocess.check_call(['sh', 'install.sh'])
        install.run(self)

setup(
    name='pylsodes',
    version='0.1.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
    ],
    cmdclass={
        'install': CustomInstallCommand,
    },
    author='S. Maitrey',
    author_email='km.maitrey@gmail.com',
    description='Python wrapper for DLSODES solver from Fortran ODEPACK family of differential equation solvers',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/kmaitreys/pylsodes',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
