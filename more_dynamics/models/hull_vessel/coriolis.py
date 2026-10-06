import casadi as ca

from .mass_properties import HullMassProperties


def _skew(vector: ca.SX) -> ca.SX:
    """Return the skew matrix whose product is ``vector x other``."""
    x, y, z = vector[0], vector[1], vector[2]
    return ca.vertcat(
        ca.horzcat(0.0, -z, y),
        ca.horzcat(z, 0.0, -x),
        ca.horzcat(-y, x, 0.0),
    )


def coriolis_matrices_casadi(
    mass_properties: HullMassProperties,
    velocity: ca.SX,
) -> tuple[ca.SX, ca.SX]:
    """Build Fossen rigid-body and added-mass Coriolis matrices.

    This is the base ``C_RB + C_A`` decomposition used by the modular
    ``HullUSV`` reference model before its optional added-mass suppression.
    Both matrices are derived from the single configured rigid-body/added-mass
    pair.
    """
    linear_velocity = velocity[:3]
    angular_velocity = velocity[3:]
    center_of_gravity = ca.DM(mass_properties.center_of_gravity)
    transform = ca.vertcat(
        ca.horzcat(ca.DM.eye(3), _skew(center_of_gravity).T),
        ca.horzcat(ca.DM.zeros(3, 3), ca.DM.eye(3)),
    )
    inertia_at_cg = ca.DM(mass_properties.inertia_at_center_of_gravity)
    rigid_body_at_cg = ca.vertcat(
        ca.horzcat(
            mass_properties.mass * _skew(angular_velocity),
            ca.DM.zeros(3, 3),
        ),
        ca.horzcat(
            ca.DM.zeros(3, 3),
            -_skew(inertia_at_cg @ angular_velocity),
        ),
    )
    rigid_body_coriolis = transform.T @ rigid_body_at_cg @ transform

    added_mass = ca.DM(mass_properties.added_mass_matrix)
    added_linear_momentum = (
        added_mass[:3, :3] @ linear_velocity
        + added_mass[:3, 3:] @ angular_velocity
    )
    added_angular_momentum = (
        added_mass[3:, :3] @ linear_velocity
        + added_mass[3:, 3:] @ angular_velocity
    )
    linear_momentum_skew = _skew(added_linear_momentum)
    angular_momentum_skew = _skew(added_angular_momentum)
    added_mass_coriolis = ca.vertcat(
        ca.horzcat(ca.DM.zeros(3, 3), -linear_momentum_skew),
        ca.horzcat(-linear_momentum_skew, -angular_momentum_skew),
    )
    return rigid_body_coriolis, added_mass_coriolis
