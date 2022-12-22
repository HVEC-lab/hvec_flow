"""
Package hydraulic engineers dealing with flow and turbulence.

Created by Hessel Voortman Engineering Consultancy, 2022
"""
from .admin import __author__, __author_email__, __version__

"""
The constants and functions in the following subpackage 
are so frequently used that they are made available directly
in the package
"""
from .constants import *
from .channelflow import *
from .energyloss import *
from .geometry import *
from .hydrostatics import *

"""
The following two remain in their own namespace
"""
from . import sills
from . import culverts

"""
The subpackage turbulence is highly specialised and
only imported when called for.
"""