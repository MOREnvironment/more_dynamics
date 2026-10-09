"""The choices every vehicle type offers beside its own parameters: how the
ocean current is taken into the body axes (``current_form``) and where the
water density and the acceleration of gravity come from (``site_form``).
One list per choice, so the three vehicle types name the same forms.

``current_form`` (fidelity, what raises it):

* ``"full_attitude"`` (default): the current constant in NED taken into the
  body axes with the full attitude and its rate (Fossen 2011, eqs.
  8.138-8.141, 8.157, pp. 221-225). Raise: a measured current profile.
* ``"none"``: still water.
* ``"yaw_rate_terms"``, ``"full_rotation_rate"``: the two shortcuts of MSS
  (``remus100.m`` 118-122; ``otter.m`` 113-118), kept to reproduce MSS; each
  is exact only for a level craft (``models/current``).

``site_form``:

* ``"latitude"``: gravity from WGS-84 normal gravity at the vehicle's
  ``latitude`` (``remus100.m`` 96-97; ``gravity.m`` 11-12).
* ``"given"``: ``gravity`` as a value, 9.81 m/s^2 to match MSS from release
  2.0.4 (``LIBRARY/mssConstants.m`` 14 @ 49e03e3; ``otter.m`` 90-91).

The water density and the kinematic viscosity are parameters of the vehicle
in both forms.

References
----------
[Fossen 2011] Fossen, T. I. (2011). Handbook of Marine Craft Hydrodynamics
    and Motion Control. Wiley. Eqs. 8.138-8.141, 8.157, pp. 221-225.
[MSS] Fossen, T. I. MSS, MIT: CRAFT/AUV/models/remus100.m 96-97, 118-122 and
    CRAFT/USV/models/otter.m 90-91, 113-118 @ cc07579; LIBRARY/mssConstants.m
    12-14 @ 49e03e3 (release 2.0.5).

Author:    Enio Krizman
Date:      2026-10-09
"""
from more_dynamics.models.current.current import (
    horizontal_current_full_rotation_rate, horizontal_current_yaw_rate_terms, no_current,
    uniform_current_full_attitude)
from more_dynamics.models.site.site import (
    SITE_AT_LATITUDE_PARAMETERS, SITE_GIVEN_PARAMETERS, site_at_latitude, site_given)

CURRENT_FORMS = {"full_attitude": uniform_current_full_attitude, "none": no_current,
                 "yaw_rate_terms": horizontal_current_yaw_rate_terms,
                 "full_rotation_rate": horizontal_current_full_rotation_rate}
CURRENT_INPUTS = (("current_speed", 1, "horizontal current speed, m/s"),
                  ("current_direction", 1, "current set, NED, from north, clockwise, rad"),
                  ("current_vertical_speed", 1, "vertical current speed, NED, down positive, m/s"))
SITE_FORMS = {"latitude": (site_at_latitude, SITE_AT_LATITUDE_PARAMETERS),
              "given": (site_given, SITE_GIVEN_PARAMETERS)}


def current_inputs(context, plugin, error):
    """The vehicle's own inputs for the chosen ``current_form`` (none for still water); ``error`` is raised for an
    unknown form."""
    form = context.get_parameter("current_form")
    if form not in CURRENT_FORMS:
        raise error(f"{plugin}: current_form must be one of {sorted(CURRENT_FORMS)}, got {form!r}")
    return () if form == "none" else CURRENT_INPUTS


def select_site(parts):
    """``(function, declared parameters)`` of the site form named by the vehicle's ``site_form`` (``Parts.select``
    refuses an unknown one, naming the plugin)."""
    return parts.select("site_form", SITE_FORMS)
