#!/usr/bin/env python

from setuptools import find_packages, setup

kwds = {'name': 'asap',
        'version': '0.0.1',
        'packages': find_packages('.'),
        'package_dir': {'asap': 'asap'},
        'author': 'Laughlin Research, LLC',
        'author_email': 'trevor@laughlinresearch.com',
        'maintainer': 'Laughlin Research, LLC',
        'maintainer_email': 'trevor@laughlinresearch.com',
        'description': 'Aerospace Structural Analysis Program',
        'license': 'Proprietary'}

setup(**kwds)
