"""A block as a member of the ``force_producers`` list slot: its own command
input is renamed ``command`` (the plugin layer gives it a payload name) and
only the listed outputs are kept.

Author:    Enio Krizman
Date:      2026-10-08
"""

from more_dynamics.models.shared.wiring import function_from


def as_producer(block, command_name, name, keep=("tau",)):
    """``block`` as a ``ForceProducer`` leaf: its own command input renamed
    ``command``."""
    ins = {("command" if n == command_name else n): block.sx_in(i) for i, n in enumerate(block.name_in())}
    out = block.call({n: ins["command" if n == command_name else n] for n in block.name_in()})
    return function_from(name, ins, {k: out[k] for k in keep})
