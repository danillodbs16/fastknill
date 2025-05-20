from setuptools import setup, Extension

# Define your Cython extension modules here
extensions = [
    Extension("fastknill", ["fastknill.pyx"]),
    # Add more extensions if you have multiple Cython modules
]

setup(
    name='fastknill',
    author='Danillo Barros de Souza',
    author_email='danillo.dbs16@gmail.com',
    description='fastknill - An alternative fast computation of Euler charactetistics and Knill curvature from networks',
    url='https://github.com/danillodbs16/fastknill',
    version='0.1.0',
    ext_modules=extensions,
    install_requires=[
    'numpy',
    'cython',
    'scikit-learn',
    'scipy'
    ],
)
