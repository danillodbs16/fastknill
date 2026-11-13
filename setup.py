
from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        name="fastknill.fastknill",
        sources=["src/fastknill/fastknill.pyx"],
    )
]

setup(
    name="fastknill",
    packages=["fastknill"],
    package_dir={"": "src"},
    ext_modules=cythonize(
        extensions,
        compiler_directives={"language_level": "3"},
    ),
)
