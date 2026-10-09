from __future__ import annotations

# REMUS 100 restoring: neutral buoyancy, W = B (remus100.m 214).
# Source of every value: the comment beside it (the plugin's own default).
# rpp rewrites this file without comments when parameters are saved from the editor.


class ComponentParameters:
    buoyancy_method = "neutral"  # remus100.m:214 B = W
    center_of_buoyancy = [0, 0, 0]  # remus100.m:138 r_bB
    displaced_volume = 0.030243065278557742  # m^3, 4/3 pi (L/2)(D/2)^2, remus100.m:131-132
    buoyancy = 313.31493718007266  # N, m g with m = 31.9 kg (remus100.m:3), g at remus100.m:96-97
