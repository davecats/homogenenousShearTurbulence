#!/usr/bin/env python3
"""
Convert the new simulation input deck (hst.in, a Fortran namelist) into the
flat "name = value" token file (readinput.in) that the CPL post-processor
(post/cpl/readinput.cpl) reads.

WHY: the new simulator (src/hst_input.f90 / src/hst_io.f90) reads a namelist
(hst.in, groups &mesh / &physics / &time_control).  The CPL post-processor's
READ BY NAME (scanrec) only accepts flat `name = value` tokens -- it halts on
the `&group` / `/` Fortran namelist delimiters, so it cannot read hst.in
verbatim.  This script is the compatibility boundary: it pulls every value the
post-processor needs out of hst.in and emits them as flat tokens.

Axis convention: the post-processor uses the convention of the simulator,
  x streamwise, y vertical (the shear direction), z spanwise,
so nx, ny, nz (modes in x, points in y, modes in z) and ly are copied from
hst.in &mesh without any relabelling.
readinput.cpl declares nx/ny/nz as plain INTEGERs (no CONSTANT), so CPL builds
its field arrays with RUNTIME extents taken from these readinput.in tokens
(CPL's compiler automatically heap-allocates arrays whose bounds are not
compile-time constants).  This script therefore emits nx/ny/nz FIRST, so
readinput.cpl reads them before touching any array.

Usage:
    python3 convert_hstin.py [hst.in] [readinput.in]
Defaults: input hst.in, output readinput.in (next to this script).
"""

import re
import sys
from pathlib import Path

# Fortran logical -> the YES/NO spelling the CPL boolean reader understands.
BOOL = {'.true.': 'YES', '.false.': 'NO', 'true': 'YES', 'false': 'NO',
        't': 'YES', 'f': 'NO'}

# Order matters: readinput.cpl issues READ BY NAME in this exact order (the
# CPL file cursor advances forward).  (name, new_deck_key, group, mapper)
TOKENS = [
    ('alfa0',            'alfa0',          'mesh'),
    ('beta0',            'beta0',          'mesh'),
    ('ly',               'ly',             'mesh'),
    ('ystretch',         'ystretch',       'mesh'),
    ('re',               're',             'physics'),
    ('S',                'S',              'physics'),
    ('deltat',           'deltat',         'time_control'),
    ('cflmax',           'cflmax',         'time_control'),
    ('t_max',            't_max',          'time_control'),
    ('dt_stat',          'dt_stat',        'time_control'),
    ('dt_field',         'dt_field',       'time_control'),
    ('dt_save',          'dt_save',        'time_control'),
    ('time_from_restart','time_from_restart','time_control'),
    ('restart_file',     'restart_file',   'time_control'),
]

def parse_deck(path):
    """Return dict key -> value_string for the given namelist deck.

    Handles multiple `name = value` pairs on one line (e.g. `nx = 65, ny = 86,
    nz = 22`) by splitting each line on commas.  Tolerates a trailing comma on
    the last pair and a trailing `! comment`.  Note: this assumes no rhs value
    itself contains a comma (true for this deck).  Quoted string values are
    kept verbatim, so readinput.cpl's STRING reader sees the unquoted name from
    convert_hstin.py after quoting is stripped there.
    """
    vals = {}
    group = None
    for line in path.read_text().splitlines():
        stripped = line.split('!', 1)[0].strip()   # drop trailing comments
        if not stripped:
            continue
        low = stripped.lower()
        if low.startswith('&'):
            group = low[1:].strip().rstrip('&').strip()
            continue
        if stripped == '/' or (low == '&end' and group):
            group = None
            continue
        for piece in stripped.split(','):
            m = re.match(r'([A-Za-z_]\w*)\s*=\s*(.*)', piece.strip())
            if m:
                vals[m.group(1)] = m.group(2).strip()
    return vals

def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('hst.in')
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).with_name('readinput.in')

    vals = parse_deck(src)

    out = []
    # Array extents come from &mesh and are emitted FIRST, verbatim: modes in x,
    # points in y (shear direction), modes in z (spanwise).  readinput.cpl reads
    # these into its plain INTEGER nx/ny/nz before any array is built, so CPL
    # sizes the field arrays at runtime.
    for key in ('nx', 'ny', 'nz'):
        out.append(f'{key} = {int(vals[key])}')

    for name, key, grp in TOKENS:
        if key not in vals:
            print(f'WARNING: {src.name} has no `{key}` in &{grp}; leaving it unset.', file=sys.stderr)
            continue
        val = vals[key]
        low = val.lower()
        if low in BOOL:
            val = BOOL[low]
        elif val.startswith("'") and val.endswith("'"):
            val = val[1:-1]                       # strip quotes for CPL string read
        out.append(f'{name} = {val}')

    dst.write_text('\n'.join(out) + '\n')
    print(f'Wrote {dst} ({len(out)} tokens) from {src}.')

if __name__ == '__main__':
    main()
