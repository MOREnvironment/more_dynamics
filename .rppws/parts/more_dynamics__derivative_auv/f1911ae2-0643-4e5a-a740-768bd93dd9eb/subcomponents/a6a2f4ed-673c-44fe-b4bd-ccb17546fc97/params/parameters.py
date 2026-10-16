from __future__ import annotations


class ComponentParameters:
    buoyancy_method = 'neutral'  # neutral: B = W as npsauv.m 214 gives it
    center_of_buoyancy = [0, 0, 0]  # npsauv.m:107 [xB yB zB] (via MSS, Healey & Lienard 1993 not read)
    displaced_volume = []  # unset: read by no selected form of this part
    buoyancy = 53400  # npsauv.m:108 B = 53400 (neutral, W = B as given) (via MSS, Healey & Lienard 1993 not read)
