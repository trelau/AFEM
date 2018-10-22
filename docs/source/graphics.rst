Graphics
========
The ``afem.graphics`` package provides a base class for viewable objects as
well as a minimal tool to visualize AFEM shapes, geometry, and meshes. The most
used tool will be the :class:`.Viewer` class which can be imported by::

    from afem.graphics import Viewer

An instance can be created, entities added, and then viewed similar to the
following process:

.. code-block:: python

    gui = Viewer()
    gui.add(*args)
    gui.start()

The ``add()`` method will try and process the argument given its type and is the
most generic method to call. More specific methods giving the user more control
are described in the class documentation. When the ``start()`` method is called
an application is launched to view the contents and program execution will stop.

The viewer instance will be destroyed if exited using the "exit" button in the
top right. If the user wants to continue processing with the same instance,
simply press the ``c`` key on the keyboard to continue processing. Hotkeys are
available as shown:

.. table:: Hotkeys for ``Viewer``
   :widths: auto

   ===== =======================================================================
   Key   Description
   ===== =======================================================================
   ``c`` Continue processing.
   ``f`` Fit the contents.
   ``s`` Shaded view.
   ``w`` Wireframe view.
   ``i`` Isometric view.
   ``t`` Top view.
   ===== =======================================================================

.. automodule:: afem.graphics.display
