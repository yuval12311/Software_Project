from setuptools import setup, find_packages, Extension

setup (
    name='kmeans',
    version='0.1',
    packages=find_packages(),
    ext_modules= [
        Extension(
            'kmeans',
            ['kmeans.c']
        )
    ]
)