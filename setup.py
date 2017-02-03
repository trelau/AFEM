#!/usr/bin/env python

import setuptools
from numpy.distutils.core import setup


def setup_asap():
    kwds = {'name': 'ASAP',
            'version': '0.0.1',
            'package_dir': {'ASAP': 'asap'},
            'packages': setuptools.find_packages('.'),
            'author': 'Laughlin Research, LLC',
            'author_email': 'trevor@laughlinresearch.com',
            'maintainer': 'Laughlin Research, LLC',
            'maintainer_email': 'trevor@laughlinresearch.com',
            'description': 'Aerospace Structural Analysis Program',
            'license': 'Proprietary'}

    setup(**kwds)


if __name__ == '__main__':
    setup_asap()
