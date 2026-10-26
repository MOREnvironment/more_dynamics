from __future__ import annotations

# Grethe (NTNU Mariner 5 USV, Maritime Robotics build 40402, Pioner 17 ft hull): hull, mass
# and forms. Source of every value: the comment beside it, reproduced from
# tests/vehicle_models/data/grethe/grethe_parameters.json (the one hand-edited source, R10,
# E-104) -- that file and this part's own SOURCE.md carry the full citation and kind.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    length = 5.2  # document: A-43 M1, F003 p2, M001 p7 (overall length)
    beam = 2.15  # document: A-43 M2, F003 p2, M001 p7, D002 p1
    hull_block_coefficient = 0.234  # derived: A-43 M4, 806 / (1025 * 5.2 * 2.15 * 0.3), a mass knob not hull fullness
    hull_waterplane_coefficient = 0.41  # fitted: grethe.m 242 Cw, chosen so that the roll period is the logged 1.85 s (a hard-chine hull would be 0.7-0.8; it absorbs the unexplained roll stiffness)
    hull_mass = 806.0  # document: A-43 M7, F003 p2 (FAT, +/- 20 kg)
    payload_mass = 0.0  # a bare hull; the 2026 field payload (USBL mount, battery packs) is not weighed (A-43 M8, gap)
    hull_center_of_gravity = [0.0, 0.0, 0.025]  # estimate: A-43 M9
    payload_position = [0.0, 0.0, 0.0]  # unused while payload_mass is 0 (A-43 M10)
    radii_of_gyration = [0.35, 0.25, 0.25]  # estimate: A-43 M12
    added_mass_coefficients = [-1.0, -1.5, -1.0, -0.2, -0.8, -1.7]  # estimate: A-43 M11 (open: library yaw is -1.2)
    added_mass_form = "scaled_derivatives"  # A-43 M11
    added_mass_mass_basis = "displaced"  # owner decision: E-101 Q5 a
    coriolis_form = "kirchhoff_full"  # every added-mass Coriolis term kept
    current_form = "none"  # this part is gated in still water; the G1 references hold current fixed
