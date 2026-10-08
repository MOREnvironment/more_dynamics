PropellerLinearizedOpenWater: MSS form | linearised open-water characteristic with given coefficients (remus100.m 148-176) | assumes the printed coefficients | raise: a measured open-water curve or the Wageningen series

value | kind | place
propeller_diameter | 0.14 | published | remus100.m:148 D_prop
max_shaft_speed | 1525 | published | remus100.m:110 n_max
thrust_deduction | 0.1 | published | remus100.m:149 t_prop
wake_fraction | 0.05600000000000005 | published | remus100.m:150 Va = 0.944 * U_r, w = 1 - 0.944
pitch_diameter_ratio | 1 | published | remus100.m:156 wageningen(0,1,0.718,3), second argument
blade_area_ratio | 0.718 | published | remus100.m:155-156 blade-area ratio 0.718
blade_count | 3 | published | remus100.m:155-156 3 blades
max_advance_number | 0.6632 | published | remus100.m:153 Ja_max
roll_moment_scale | 0.1 | published | remus100.m:252 tau(4) = K_prop / 10
position | [0, 0, 0] | published | remus100.m:249-252 thrust on x_b through the CO (no moment arm)
orientation | [0, 0, 0] | published | remus100.m:249, 252 shaft along x_b (thrust in tau(1), torque in tau(4))
thrust_torque_coefficients | [0.4566, 0.07, 0.1798, 0.0312] | published | remus100.m:157, 158, 160, 161 [KT_0 KQ_0 KT_max KQ_max]
