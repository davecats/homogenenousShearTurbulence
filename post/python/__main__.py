"""Command-line dispatcher for the postprocessing tools.

Usage::

    python3 -m post.python <tool> [deck] [tool args...]

``deck`` is the run directory or an ``hst.in`` path (default: ``hst.in`` in the
current working directory).  The remaining arguments are passed to the tool.
List the available tools with ``python3 -m post.python --list``.
"""

from __future__ import annotations

import sys


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if argv and argv[0] in ("--list", "-l", "list"):
        from .tools import TOOLS
        print("available postprocessing tools:")
        for name in sorted(TOOLS):
            print(f"  {name}")
            doc = (TOOLS[name].__doc__ or "").strip().splitlines()
            if doc:
                print(f"      {doc[0]}")
        return 0

    if not argv:
        print(__doc__)
        return 1

    tool_name = argv.pop(0)
    if not tool_name.startswith("-"):
        # the first non-flag token may be a tool; the deck is then the next
        pass

    # The tool name is the first argument.  The deck is either the next
    # positional that looks like a path/directory, or defaults to hst.in.
    # We pass everything after the tool name to the tool, which decides how
    # much is its own; the deck is taken as the first positional that is not
    # consumed as a flag value.
    from .tools import TOOLS
    tool = TOOLS.get(tool_name)
    if tool is None:
        print(f"unknown tool '{tool_name}'", file=sys.stderr)
        print("run `python3 -m post.python --list` for the available tools", file=sys.stderr)
        return 2

    # Determine the deck.  An explicit first positional that looks like a run
    # directory or a namelist file is consumed and used verbatim.  Otherwise we
    # look for ``hst.in`` in the current directory, then walk up from each
    # non-flag positional (so ``slices run/fields/field1.fld ...`` finds
    # ``run/hst.in`` for the headerless pressure files).
    from pathlib import Path
    from .config import load_config

    deck = None
    if argv and not argv[0].startswith("-"):
        cand = Path(argv[0])
        if cand.is_dir() or cand.name == "hst.in" or (cand.is_file() and cand.suffix == ".in"):
            deck = argv.pop(0)

    if deck is not None:
        cfg = load_config(deck)
        return tool.run(cfg, argv)

    # No explicit deck: prefer ``hst.in`` in the current directory, else the
    # deck of the (first) run directory that a positional path lives in.
    cfg = None
    try:
        cfg = load_config(None)          # hst.in in CWD
    except FileNotFoundError:
        cfg = None
    if cfg is None and argv:
        for token in argv:
            if token.startswith("-"):
                continue
            base = Path(token)
            if not base.exists():
                continue
            p = base if base.is_dir() else base.parent
            for _ in range(6):
                if (p / "hst.in").is_file():
                    cfg = load_config(p / "hst.in")
                    break
                if p == p.parent:
                    break
                p = p.parent
            if cfg is not None:
                break
    if cfg is None:
        cfg = load_config(None)          # raise the canonical FileNotFoundError

    return tool.run(cfg, argv)


if __name__ == "__main__":
    sys.exit(main())
