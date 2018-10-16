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
from afem.smesh.entities import MeshGen
from afem.smesh.hypotheses import (Regular1D, NetgenAlgo2D,
                                   NetgenSimple2D, LocalLength1D)
from afem.structure.group import GroupAPI

__all__ = ["MeshVehicle"]


class MeshVehicle(object):
    """
    Tool for assisting in vehicle-level mesh generation. The "master" group
    will be used to define the top-level shape for meshing, which should
    include all the parts and therefore shapes to be meshed. Default mesh
    controls are applied to the top-level shape based on the target element
    size.

    :param float target_size: Default global element size.
    :param bool allow_quads: Option to generate quad-dominated mesh.
    """

    def __init__(self, target_size, allow_quads=True):
        group = GroupAPI.get_master()
        self._shape = group.prepare_shape_to_mesh()
        self._gen = MeshGen()
        self._mesh = self._gen.create_mesh(self._shape)

        # Initialize each part for meshing
        for part in group.get_parts():
            part.init_meshing(self._mesh)

        # Define global mesh control based on target size
        hyp1d = LocalLength1D(self._gen, target_size)
        alg1d = Regular1D(self._gen)
        hyp2d = NetgenSimple2D(self._gen, target_size, allow_quads=allow_quads)
        alg2d = NetgenAlgo2D(self._gen)
        self.add_controls([hyp1d, alg1d, hyp2d, alg2d])

    @property
    def shape(self):
        """
        :return: The top-level shape.
        :rtype: afem.topology.entities.Shape
        """
        return self._shape

    @property
    def gen(self):
        """
        :return: The mesh generator.
        :rtype: afem.smesh.entities.MeshGen
        """
        return self._gen

    @property
    def mesh(self):
        """
        :return: The top-level mesh.
        :rtype: afem.smesh.entities.Mesh
        """
        return self._mesh

    def add_control(self, control, shape=None):
        """
        Add a mesh control.

        :param afem.smesh.hypotheses.Hypothesis control: The control.
        :param shape: The shape the control applies to. This can be a
            sub-shape of the master shape. If not provided then the master
            shape is used.
        :type shape: afem.topology.entities.Shape

        :return: Status of adding hypothesis.
        :rtype: OCCT.SMESH.SMESH_Hypothesis.Hypothesis_Status
        """
        if shape is None:
            shape = self.shape

        return self._mesh.add_hypothesis(control, shape)

    def add_controls(self, controls, shape=None):
        """
        Add mesh controls to the shape.

        :param controls: The controls.
        :type controls:
            collections.Sequence(afem.smesh.hypotheses.Hypothesis)
        :param shape: The shape the control applies to. This can be a
            sub-shape of the master shape. If not provided then the master
            shape is used.
        :type shape: afem.topology.entities.Shape

        :return: Dictionary where the key is the control and the value is
            the control status.
        :rtype: dict

        :raise ValueError: If no shape is available to apply the control to.
        """
        status_dict = {}
        for hyp in controls:
            status = self.add_control(hyp, shape)
            status_dict[hyp] = status
        return status_dict

    def compute(self):
        """
        Compute the mesh.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        return self._gen.compute(self._mesh, self._shape)
