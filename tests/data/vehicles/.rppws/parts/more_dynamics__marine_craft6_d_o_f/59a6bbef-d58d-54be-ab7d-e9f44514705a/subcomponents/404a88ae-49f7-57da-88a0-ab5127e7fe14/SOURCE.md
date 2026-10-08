HullLiftDrag: physics form | slender-body lift and induced drag (forceLiftDrag.m) | assumes the printed planform and Oswald factor | raise: towing-tank or log identification

value | kind | place
planform_area | 0.2128 | derived | remus100.m:133 S = 0.7 * L_auv * D_auv
parasitic_drag_coefficient | 0.05595961914206819 | derived | remus100.m:143-144 CD_0 = Cd * pi * b^2 / S, b = D_auv / 2 (one geometry)
oswald_efficiency | 0.3 | published | coeffLiftDrag.m:54 e = 0.3
