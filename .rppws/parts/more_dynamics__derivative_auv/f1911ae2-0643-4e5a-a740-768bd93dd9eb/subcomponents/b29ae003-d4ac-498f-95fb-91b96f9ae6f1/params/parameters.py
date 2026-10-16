from __future__ import annotations


class ComponentParameters:
    site_form = 'given'  # the given form: gravity and density as npsauv.m states them
    latitude = []  # unset: read by no selected form of this part
    gravity = 9.81  # npsauv.m:105 g = 9.81 (via MSS, Healey & Lienard 1993 not read)
    water_density = 1025  # npsauv.m:109 rho = 1025 (via MSS, Healey & Lienard 1993 not read)
    kinematic_viscosity = 1e-06  # NOT AN NPS VALUE and read by no NPS form: npsauv.m has no viscosity; the Site plugin declares it, so a build needs one; 1e-06 m^2/s is the value cylinderDrag.m 78-80 and the REMUS 100 / Otter parts use; owner to approve or replace
