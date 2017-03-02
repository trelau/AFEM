# # Add local Netgen to path if available.
# import os
# import sys
#
#
# def _add_netgen_path():
#     # Only win64 binaries are available for now.
#     if sys.platform != 'win32':
#         return False
#
#     fdir = os.path.dirname(__file__)
#     if sys.version_info[:2] == (2, 7):
#         netgen_path = '/netgen/py27_win64_msvc2015/lib/site-packages'
#     elif sys.version_info[:2] == (3, 5):
#         netgen_path = '/netgen/py35_win64_msvc2015/lib/site-packages'
#     else:
#         return False
#     path_name = '/'.join([fdir, netgen_path])
#     sys.path.append(path_name)
#     try:
#         # noinspection PyUnresolvedReferences
#         import netgen
#         return True
#     except ImportError:
#         return False
#
#
# try:
#     # noinspection PyUnresolvedReferences
#     import netgen
# except ImportError:
#     status = _add_netgen_path()
