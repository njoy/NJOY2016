from importlib.metadata import version, PackageNotFoundError
from .paths import *

try:
    __version__ = version('njoy')
except PackageNotFoundError:
    __version__ = "unknown"