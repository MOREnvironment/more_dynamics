from __future__ import annotations

# REMUS 100 (estimated) restoring: neutral buoyancy, W = B.
# Source of every value: the comment beside it.
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    buoyancy_method = "neutral"  # B = W
    center_of_buoyancy = [0, 0, 0]  # m, CO -> CB
    displaced_volume = 0.030243065278557742  # m^3, 4/3 pi (L/2)(D/2)^2 (derived)
    buoyancy = 313.31493718007266  # N, m g with m = 31.9 kg (derived)
