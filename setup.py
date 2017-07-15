#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
        name='asap',
        version='0.1.0',
        packages=find_packages('.'),
        package_dir={'asap': 'asap'},
        author='Laughlin Research, LLC',
        author_email='info@laughlinresearch.com',
        description='Aerospace Structural Analysis Program',
        license='Proprietary',
        install_requires=['numpy', 'scipy']
)
