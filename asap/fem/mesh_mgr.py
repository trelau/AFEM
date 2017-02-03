class MeshMgr(object):
    """
    Mesh data.
    """
    _shape_to_mesh = {}

    @classmethod
    def shape_to_mesh(cls, shape, mesh):
        """
        Associate the mesh to the shape.

        :param shape:
        :param mesh:

        :return:
        """
        cls._shape_to_mesh[str(shape)] = mesh
        return True

    @classmethod
    def mesh_from_shape(cls, shape):
        """
        Get the mesh associated with the shape.

        :param shape:

        :return:
        """
        try:
            return cls._shape_to_mesh[str(shape)]
        except KeyError:
            return None

    @classmethod
    def has_mesh(cls, shape):
        """
        Check to see if a shape has a mesh.

        :param shape:

        :return:
        """
        return str(shape) in cls._shape_to_mesh
