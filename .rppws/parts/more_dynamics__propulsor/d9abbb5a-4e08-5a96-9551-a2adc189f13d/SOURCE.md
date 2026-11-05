# REMUS 100 (estimated) propulsor

A named part of `more_dynamics::Propulsor`: the numbers of `REMUS 100 propulsor` (the open-water coefficients and geometry printed for the vehicle, see that part), with the shaft-torque reaction on the hull at the estimated scale 0.35 in place of 0.1. "REMUS 100 (estimated)" links it as its propulsor (the shaft speed command `drive_command`, rpm).

`shaft_torque_scale` 0.35: kind **estimate**. The reaction on the hull is the shaft torque times this scale (Prestero 2001, Sec. 4.6.2, p. 33); the value was calibrated in the generic-models REMUS actuator set to a steady roll of -5.3 deg at 1500 rpm and 1.51 m/s with that set's plant and is carried over here, not re-calibrated for this plant. Every other value is that of `REMUS 100 propulsor`.
