# This file is part of AFEM which provides an engineering toolkit for airframe
# finite element modeling during conceptual design.
#
# Copyright (C) 2016-2018  Laughlin Research, LLC (info@laughlinresearch.com)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
from setuptools import setup

setup(
    name='afem',
    version='0.3.0',
    packages=['afem', 'afem.adaptor', 'afem.base', 'afem.exchange', 'afem.fem',
              'afem.geometry', 'afem.graphics', 'afem.misc', 'afem.occ',
              'afem.oml', 'afem.sketch', 'afem.smesh', 'afem.structure',
              'afem.topology'],
    author='Laughlin Research, LLC',
    author_email='info@laughlinresearch.com',
    description='Airframe Finite Element Modeler',
    license='LGPL v2.1'
)
