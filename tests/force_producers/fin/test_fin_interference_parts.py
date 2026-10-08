"""Gate tests of the fin interference parts (hull-fin interference ->
factors on the deflection and on the flow angle): ``none`` now,
``slender_body`` written and waiting.

Written 2026-10-08, before the parts exist (contract: ``fin_parts_contract``).
No coupling input; outputs ``deflection_factor`` and ``flow_angle_factor``
(1 each, dimensionless). The skeleton forms ``alpha = deflection_factor *
deflection - flow_angle_factor * flow_angle``.

* ``none``: both factors 1 (MSS ``remus100.m`` 238-245 and Prestero 2001
  eqs. 4.41-4.43, pp. 32-33, apply no factor).
* ``slender_body``: factors from body radius over fin semispan (Pitts 1957,
  through a citing memo; **unread** by this job). Not gated until the source
  is read in full; the test that waits for it is below.

Prestero's hull effect on the fin is not a factor on the angle: it doubles
the aspect ratio (``AR_e = 2 AR``, eq. 4.39, p. 31) inside the lift slope,
so it belongs to the section part ``lifting_line``
(``test_fin_section_parts``).

References
----------
[MSS] Fossen, T. I. (2026). *Marine Systems Simulator (MSS)*, release 2.0.2
    with the fixes of 2026-10-07. https://github.com/cybergalactic/MSS, MIT
    licence, revision ``cc07579``: ``CRAFT/AUV/models/remus100.m`` 238-245.
[Prestero 2001] Prestero, T. (2001). *Verification of a six-degree of
    freedom simulation model for the REMUS autonomous underwater vehicle*.
    MSc thesis, MIT/WHOI. Ch. 4, eqs. 4.39, 4.41-4.43, pp. 31-33.

Author:    Enio Krizman
Date:      2026-10-08
"""

import pytest

from fin_parts_contract import WAITING_REASON, call, part, part_declared


def test_no_interference_declares_no_parameter_and_has_no_input():
    assert tuple(part_declared("interference", "none")) == ()
    assert part("interference", "none").name_in() == []


def test_no_interference_factors_are_one():
    out = call(part("interference", "none"))
    assert out["deflection_factor"][0] == 1.0 and out["flow_angle_factor"][0] == 1.0


@pytest.mark.skip(reason=WAITING_REASON)
def test_slender_body_factors_tend_to_one_as_the_body_radius_vanishes():
    """Slender-body interference -> ``none`` when body radius / fin semispan
    -> 0 (the reduction arrow of the design note's part-swap table). Written
    now; its parameter names and form are fixed when Pitts 1957 is read."""
    raise AssertionError("contract not fixed")
