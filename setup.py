#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
        name='afem',
        version='0.1.0',
        packages=find_packages('.'),
        package_dir={'afem': 'afem'},
        author='Laughlin Research, LLC',
        author_email='info@laughlinresearch.com',
        description='Aiframe Finite Element Modeling',
        license='Proprietary',
        install_requires=['numpy', 'scipy']
)
