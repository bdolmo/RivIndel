from setuptools import setup, find_packages
from Cython.Build import cythonize


with open("README.md") as fh:
    long_description = fh.read()

setup(
    name="rivindel",
    spython_requires=">=3.6",
    version="0.0.1",
    author="Bernat del Olmo",
    author_email="bernatdelolmo@gmail.com",
    description="RivIndel: Complex indel detection from short read sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="complex indel detection",
    url="https://github.com/bdolmo/RivIndel",
    install_requires=[
        "edlib>=1.3.9",
        "requests>=2.18.4",
        "pysam>=0.16.0.1",
        "Cython>=0.29.23",
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    package_dir={'': "src"},
    packages=["rivindel"],
    ext_modules=cythonize(["src/rivindel/assembler.pyx"]),

)
