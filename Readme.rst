Package hvec_flow. Contains Python implementation of
frequently used formulae for flow in hydraulic engineering.

Subpackages:

- channelflow; functions related to open channel flow
- constants; list of regularly used constants
- culverts; functions for flow in culverts
- energyloss; functions for calculating flow friction and loss coefficients
- geometry; functions describing flow geometry
- hydrostatic; functions for hydrostatic properties
- sills; functions for flow over sills
- turbulence; functions for calculation of turbulence

Prepared by Hessel Voortman Engineering Consultancy, 2022.


Installation
------------

.. code-block:: console

    pip install git+https://github.com/HVEC-lab/hvec_flow

From a clone of the repository: ``uv sync``.

Tests
-----

.. code-block:: console

    uv run pytest
