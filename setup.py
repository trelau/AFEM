#!/usr/bin/env python
from setuptools import setup

setup(
    name='afem',
    version='0.2.0',
    packages=['afem', 'afem.exchange', 'afem.geometry', 'afem.graphics',
              'afem.misc', 'afem.occ', 'afem.oml', 'afem.smesh',
              'afem.structure', 'afem.topology'],
    package_data={'afem': ['graphics/resources/*']},
    author='Laughlin Research, LLC',
    author_email='info@laughlinresearch.com',
    description='Airframe Finite Element Modeling',
    license='Proprietary'
)
