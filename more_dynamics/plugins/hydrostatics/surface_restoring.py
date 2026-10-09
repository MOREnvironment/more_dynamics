"""Restoring of a surface craft of one or two hulls, as Luka's
``HydrostaticsModel``: the displaced volume is the mass over the density, the
draft follows from the block coefficient of the hull, and the restoring matrix
is the metacentric one (Fossen 2011, ch. 4, eq. 4.25, p. 65; MSS ``otter.m``
121-122, 172-193). The first output is the signed generalized force the
vehicle adds (``vehicle_model.capnp`` 11-13), so ``restoring_force = -g``.

The hull geometry (``length``, ``hull_beam``, ``hull_block_coefficient``,
``hull_waterplane_coefficient``, ``hull_lateral_offset`` for two hulls) and
the mass, density, gravity and centre of gravity are the vehicle's
quantities, read by name; the plugin declares only what is its own. Outputs
beside the force: ``displaced_volume``, ``draft``, ``wetted_surface``, ``restoring_matrix`` (``G``
at the origin of the point P) and ``restoring_matrix_at_flotation``.

``wetted_surface`` (the wetted area of the hulls at that draft, read by the
hull loads' ITTC surge resistance) follows ``wetted_surface_method``:
``"computed"`` is the Mumford approximation of the hulls (MSS ``XuuITTC.m``
38), ``"given"`` takes the parameter ``wetted_surface``, ``"regression_table"``
is the published Radojcic et al. (2014) hydrostatic (zero-speed) wetted-surface
regression (one hull only; ``models/restoring/restoring_parts.py``). Raise:
the wetted surface from the hull lines.

Defaults are the Otter (MSS ``otter.m`` 177-179, 192): ``hull_count`` 2,
``longitudinal_inertia_factor`` 0.8, ``longitudinal_center_of_flotation``
-0.2 m. Fidelity: the linear metacentric restoring around equilibrium;
raise: a measured draft and the waterplane from the hull lines.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eq. 4.25, p. 65 (through the surface block,
    which cites ch. 4 line by line).
[MSS] Fossen, T. I. MSS, MIT, CRAFT/USV/models/otter.m 121-122, 172-193 and
    LIBRARY/modeling/XuuITTC.m 38 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-09
"""
from rpp_plugin_types.more_dynamics import HydrostaticsModel
from rpp_py.context import ComponentContext
from rpp_py.parameter_description import ParameterDescription

from more_dynamics.models.restoring.restoring_parts import (
    SURFACE_RESTORING_PARAMETERS, WETTED_SURFACE_PARAMETER, surface_restoring)
from more_dynamics.plugins.shared.payload_io import PayloadBuilder, frozen_block, payload_name


WETTED_SURFACE_METHODS = ("computed", "given", "regression_table")


class SurfaceRestoring(HydrostaticsModel):
    PARAMETERS = [
        ParameterDescription("hull_count", 2),  # otter.m 172-193: two pontoons
        ParameterDescription("longitudinal_inertia_factor", 0.8),  # otter.m 179 (I_L = 0.8 ...)
        ParameterDescription("longitudinal_center_of_flotation", -0.2),  # otter.m 177, 192
        ParameterDescription("reference_point", [0.0, 0.0, 0.0]),  # the body-frame origin (CO)
        ParameterDescription("wetted_surface_method", "computed"),  # S = n 1.025 L (C_b B + 1.7 T), XuuITTC.m:38 (Mumford)
        ParameterDescription("wetted_surface", 1.77),  # m^2, given value, read with "given": the Mumford value at the Otter set (otter.m 92-107, mass 80 kg, XuuITTC.m:38)
    ]

    def __init__(self) -> None:
        self._model = None

    def initialize(self, context: ComponentContext) -> None:
        hull_count = context.get_parameter("hull_count")
        if hull_count not in (1, 2):
            raise ValueError(f"hull_count must be 1 or 2, got {hull_count!r}")
        method = context.get_parameter("wetted_surface_method")
        if method not in WETTED_SURFACE_METHODS:
            raise ValueError(f"wetted_surface_method must be one of {WETTED_SURFACE_METHODS}, got {method!r}")
        declared = SURFACE_RESTORING_PARAMETERS + ((WETTED_SURFACE_PARAMETER,) if method == "given" else ())
        self._model = frozen_block(context, surface_restoring(hull_count, wetted_surface_method=method), declared)

    def graph(self) -> HydrostaticsModel.CasadyPayload:
        if self._model is None:
            raise RuntimeError("SurfaceRestoring must be initialized before graph()")
        io = PayloadBuilder(HydrostaticsModel.CasadyPayload())
        out = io.call(self._model)
        io.output("restoring_force", -out["g"], "signed generalized force added by the vehicle, BODY, N and N m")
        for name in self._model.name_out():
            if name != "g":
                io.output(payload_name(name), out[name])
        return io.payload()
