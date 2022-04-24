from setuptools import setup, find_packages, Extension

setup (
    name='spkmeans',
    version='0.1',
    packages=find_packages(),
    ext_modules= [
        Extension(
            'spkmeans',
            ['spkmeansmodule.c', 'spkmeans.c', 'kmeans.c', 'vector.c'],
        ),
    ]
)