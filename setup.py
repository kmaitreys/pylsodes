from setuptools import find_packages, setup

setup(
    name="pylsodes",
    version="0.2.1",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy<2.0",
    ],
    author="S. Maitrey",
    author_email="km.maitrey@gmail.com",
    description="Python wrapper for DLSODES solver from Fortran ODEPACK family of differential equation solvers",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/kmaitreys/pylsodes",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
    ],
    python_requires=">=3.10",
)
