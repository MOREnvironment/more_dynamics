"""The choices every vehicle type offers beside its own parameters: how the
ocean current is taken into the body axes (``current_form``) and, for the
types that offer it, how the attitude is carried (``kinematics_form``). One
list, so the vehicle types name the same forms. The water density, the gravity and
the kinematic viscosity are not a choice of the vehicle: they come from the
component in its ``site`` slot (``plugins/site/site.py``).

``current_form`` (fidelity, what raises it):

* ``"full_attitude"`` (default): the current constant in NED taken into the
  body axes with the full attitude and its rate (Fossen 2011, eqs.
  8.138-8.141, 8.157, pp. 221-225). Raise: a measured current profile.
* ``"none"``: still water.
* ``"yaw_rate_terms"``, ``"full_rotation_rate"``: the two shortcuts of MSS
  (``remus100.m`` 118-122; ``otter.m`` 113-118), kept to reproduce MSS; each
  is exact only for a level craft (``models/current``).

``kinematics_form`` (fidelity, what raises it):

* ``"euler"`` (default): the state holds ZYX Euler angles (12 states); singular at a pitch of +-90 deg.
* ``"unit_quaternion"``: the state holds the unit quaternion ``[eta, eps1, eps2, eps3]`` (13 states, no
  singularity; Fossen 2011, eq. 2.69, p. 30; ``remus100.m`` 225-227 @ ed2ce7b). The vehicle's output stays the 12
  values position, Euler angles, velocity, so the consumers of the pose do not change. The quaternion is kept at unit
  length after every integration step of the vehicle's own step.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 8.138-8.141, 8.157, pp. 221-225.
[MSS] Fossen, T. I. MSS, MIT: CRAFT/AUV/models/remus100.m 118-122 and
    CRAFT/USV/models/otter.m 113-118 @ cc07579.

Author:    Enio Krizman
Date:      2026-10-09
"""
from more_dynamics.models.current.current import (
    horizontal_current_full_rotation_rate, horizontal_current_yaw_rate_terms, no_current,
    uniform_current_full_attitude)

CURRENT_FORMS = {"full_attitude": uniform_current_full_attitude, "none": no_current,
                 "yaw_rate_terms": horizontal_current_yaw_rate_terms,
                 "full_rotation_rate": horizontal_current_full_rotation_rate}
KINEMATICS_FORMS = ("euler", "unit_quaternion")
CURRENT_INPUTS = (("current_speed", 1, "horizontal current speed, m/s"),
                  ("current_direction", 1, "current set, NED, from north, clockwise, rad"),
                  ("current_vertical_speed", 1, "vertical current speed, NED, down positive, m/s"))


def current_inputs(context, plugin, error):
    """The vehicle's own inputs for the chosen ``current_form`` (none for still water); ``error`` is raised for an
    unknown form."""
    form = context.get_parameter("current_form")
    if form not in CURRENT_FORMS:
        raise error(f"{plugin}: current_form must be one of {sorted(CURRENT_FORMS)}, got {form!r}")
    return () if form == "none" else CURRENT_INPUTS


def kinematics_form(context, plugin, error):
    """The vehicle's ``kinematics_form``; ``error`` is raised for an unknown one."""
    form = context.get_parameter("kinematics_form")
    if form not in KINEMATICS_FORMS:
        raise error(f"{plugin}: kinematics_form must be one of {list(KINEMATICS_FORMS)}, got {form!r}")
    return form
