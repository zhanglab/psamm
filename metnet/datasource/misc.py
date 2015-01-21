
"""Miscellaneaus data source functions"""

import re
from decimal import Decimal

from metnet.reaction import Reaction, Compound

def parse_sudensimple_reaction(s):
    """Parse a reaction string (SudenSimple)"""

    def parse_compound_list(s):
        for cpd in s.split('+'):
            if cpd == '':
                continue
            count, spec = cpd.strip().split(' ')
            spec_split = spec.split('[')
            if len(spec_split) == 2:
                # compartment
                cpdid = spec_split[0]
                comp = spec_split[1].rstrip(']')
            else:
                cpdid = spec
                comp = None
            d = Decimal(count)
            if d % 1 == 0:
                d = int(d)
            yield Compound(cpdid, compartment=comp), d

    cpd_left, cpd_right = s.split('<=>')
    left = parse_compound_list(cpd_left)
    right = parse_compound_list(cpd_right)

    return Reaction('<=>', left, right)

def parse_metnet_reaction(s):
    """Parser for the reaction format in MetNet model"""

    def parse_compound_list(s, global_comp):
        if s.strip() == '':
            return

        for cpd in s.strip().split('+'):
            cpdsplit = cpd.strip().split(') ')
            if len(cpdsplit) == 2:
                count = Decimal(cpdsplit[0].lstrip('('))
                if count % 1 == 0:
                    count = int(count)
                cpdspec = cpdsplit[1]
            else:
                count = 1
                cpdspec = cpd.strip()

            if global_comp is None:
                cpd_spec_split = cpdspec.split('[')
                if len(cpd_spec_split) == 2:
                    comp = cpd_spec_split[1].rstrip(']')
                    cpdid = cpd_spec_split[0]
                else:
                    cpdid = cpd_spec_split
                    comp = None
            else:
                cpdid = cpdspec
                comp = global_comp

            yield Compound(cpdid, compartment=comp), count

    # Split by colon for compartment information
    m = re.match(r'^\s*\[(\w+)\]\s*:\s*(.*)', s)
    if m is not None:
        global_comp = m.group(1)
        s = m.group(2)
    else:
        global_comp = None

    # Split by equation arrow
    direction = Reaction.Right
    left_right = s.split('-->')
    if len(left_right) == 1:
        direction = Reaction.Bidir
        cpd_left, cpd_right = s.split('<==>')
    else:
        cpd_left, cpd_right = left_right

    # Remove incompatible characters from compound id
    def translate(cpdid):
        return 'cpd_' + cpdid.replace(',', '__')

    # Split by +
    left = ((compound.translate(translate), count) for compound, count in parse_compound_list(cpd_left, global_comp))
    right = ((compound.translate(translate), count) for compound, count in parse_compound_list(cpd_right, global_comp))

    return Reaction(direction, left, right)
