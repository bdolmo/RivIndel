from setuptools import setup
from Cython.Build import cythonize

setup(
    name = "RivIndel",
	spython_requires='>=3.6',
    version = 1.0,
    author = "Bernat del Olmo",
    author_email = "bernatdelolmo@gmail.com",
    description = ("RivIndel NGS Pipeline"),
    keywords = "complex indel detection",
    ext_modules = cythonize("src/assembler.pyx"),
    url = "https://github.com/bdolmo/RivIndel",

    install_requires=[
        'edlib>=1.3.9',
        'requests>=2.18.4',
        'pysam>=0.16.0.1',
        'Cython>=0.29.23'
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

)
