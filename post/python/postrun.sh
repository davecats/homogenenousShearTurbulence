#!/usr/bin/env bash
# =============================================================================
# postrun.sh  --  configurable launcher for the `hst` postprocessing pipeline.
#
# Convenience wrapper around `python3 -m post.python <tool> ...` (see post/python/).
#
# Usage:
#     ./post/python/postrun.sh [OPTIONS] <command> [command args...]
#     (or:  post/python/postrun.sh ...   from the repo root)
#
# Commands:
#     stats              plot run statistics (Runtimedata + variances)
#     slice  <field>     extract a slice from a velocity/pressure snapshot
#     all                run everything (stats + slice of every snapshot found)
#     list               list the available pipeline tools
#     help               this help
#
# Global OPTIONS:
#     -r <dir>        run directory  (default: RUN_DIR below)
#     -d <deck>       path to an hst.in deck (overrides -r)
#     -o <dir>        output directory for images/.npz (default: OUTPUT_DIR)
#     -h              this help
#
# Options for the `slice` sub-command (after the snapshot path):
#     -c <u|v|w|p>    component (default: DEFAULT_COMPONENT)
#     -z <iz>         x-y plane at spanwise index iz (default: mid-box)
#     -y <j>          y-z plane at vertical row j (0..ny-1)
#     --npz           write a .npz instead of a .png
#     -o <file>       explicit output path (overrides the auto-derived name)
#
# Most settings are configurable at the top of the file (the "CONFIG" block);
# command-line options take precedence.
# =============================================================================
set -euo pipefail

# -----------------------------------------------------------------------------
# CONFIG  -- edit these as needed
# -----------------------------------------------------------------------------
# Directory this script lives in.  Since postrun.sh sits inside the `post/python`
# package directory, POST_ROOT is the package dir; REPO_ROOT (its parent's
# parent) is the repo root that must be on PYTHONPATH for `python3 -m post.python`
# to resolve `post.python`.
POST_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$POST_ROOT/../.." && pwd)"

# Run directory / name of a deck.  Leave RUN_DIR empty ("") to run from the
# current directory (auto-discover).  DECK overrides RUN_DIR when non-empty.
RUN_DIR="/home/ws/th5982/research/codes/homogenenousShearTurbulence/run"                 # e.g. "/path/to/run"  (top of run/)
DECK="/home/ws/th5982/research/codes/homogenenousShearTurbulence/run/hst.in"                    # e.g. "/path/to/run/hst.in"

# Where outputs are written.  Leave "" to use the run directory.
OUTPUT_DIR="/home/ws/th5982/research/codes/homogenenousShearTurbulence/run/postpro"

# Which snapshot index(es) to slice in `all` when not given explicitly.
# Set to a number like "3", a list like "1 3 5", or leave empty to slice all.
SLICE_INDEXES=""

# Defaults for `slice` (overridable on the command line / env).
DEFAULT_COMPONENT="u"      # u | v | w | p
DEFAULT_PLANE=""           # "" (mid-box x-y), "-z <iz>", or "-y <row>"
DEFAULT_OUT_KIND="png"     # png | npz; use "npz" for data-only output

# Separator used to visually split the output of multi-file runs.
SEP="------------------------------------------------------------------"
# =============================================================================

PYTHON_BIN="${PYTHON_BIN:-python3}"

# --- print the usage header (the comment block at the top of this file) ------
show_help() {
    awk 'NR==1 {next} /^#/ {print; next} {exit}' "$0" | sed 's/^# \{0,1\}//'
}


# --- parse global options ----------------------------------------------------
run_dir="$RUN_DIR"
deck="$DECK"
out_dir="$OUTPUT_DIR"
while getopts ":r:d:o:h" opt; do
    case "$opt" in
        r) run_dir="$OPTARG" ;;
        d) deck="$OPTARG" ;;
        o) out_dir="$OPTARG" ;;
        h) show_help; exit 0 ;;
        \?) echo "unknown option -$OPTARG (see $0 -h)" >&2; exit 2 ;;
    esac
done
shift $((OPTIND - 1))

[ $# -ge 1 ] || { echo "no command given (see $0 -h)" >&2; exit 2; }
command="$1"
shift

# --- PYTHONPATH so `python3 -m post.python` is importable -------------------
export PYTHONPATH="$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}"

# --- resolve the deck argument passed to the pipeline ------------------------
# Explicit DECK  >  RUN_DIR (as a directory)  >  nothing (auto-discover).
# Normalize trailing slashes so paths read ``run/fields/...`` not ``run//fields``.
run_dir="${run_dir%/}"
[ -n "$deck" ] && deck="${deck%/}"
[ -n "$out_dir" ] && out_dir="${out_dir%/}"
if [ -n "$deck" ]; then
    deck_arg="$deck"
elif [ -n "$run_dir" ]; then
    deck_arg="$run_dir"
else
    deck_arg=""
fi

# --- helpers ------------------------------------------------------------------
# resolve_out <default-basename> -> path under the chosen output directory
resolve_out() {
    local dir
    if [ -n "$out_dir" ]; then
        dir="$out_dir"
    elif [ -n "$run_dir" ]; then
        dir="$run_dir"
    else
        dir="$(pwd)"
    fi
    mkdir -p "$dir"
    printf '%s/%s' "$dir" "$1"
}

cmd_stats() {
    local out="$1"; shift
    echo "$SEP"
    echo "==> plot_stats"
    # shellcheck disable=SC2086
    "$PYTHON_BIN" -m post.python plot_stats $deck_arg -o "$out" "$@"
    echo "    wrote: $out"
}

cmd_slice() {
    local field="" comp plane kind=() o_override=""
    comp="$DEFAULT_COMPONENT"
    plane="$DEFAULT_PLANE"

    # Parse the slice sub-command options explicitly so the generated filename
    # always matches the requested component / plane / output kind.
    while [ $# -gt 0 ]; do
        case "$1" in
            -c) comp="$2"; shift 2 ;;
            -z) plane="-z $2"; shift 2 ;;
            -y) plane="-y $2"; shift 2 ;;
            --npz) kind="npz"; shift ;;
            -o) o_override="$2"; shift 2 ;;
            --) shift; break ;;
            -*) echo "slice: unknown option '$1'" >&2; exit 2 ;;
            *) [ -z "$field" ] && field="$1" || { echo "slice: too many positionals" >&2; exit 2; }; shift ;;
        esac
    done
    [ -n "$field" ] || { echo "slice needs a <field> (e.g. fields/field1.fld)" >&2; exit 2; }

    local name_kind="${kind:-$DEFAULT_OUT_KIND}"

    # default output name from the field filename, component & plane
    local stem spec zval
    stem="$(basename "${field%.fld}")"
    spec="$stem-c$comp"
    case "$plane" in
        -z*) zval="${plane#-z}"; spec="${spec}-z${zval// /}" ;;
        -y*) zval="${plane#-y}"; spec="${spec}-y${zval// /}" ;;
    esac
    local out
    out="${o_override:-$(resolve_out "slice-$spec.$name_kind")}"

    echo "$SEP"
    echo "==> slice  field=$field  comp=$comp  plane=[${plane:-midbox}]  kind=$name_kind"
    # shellcheck disable=SC2086
    "$PYTHON_BIN" -m post.python slices $deck_arg "$field" $plane -c "$comp" -o "$out"
    echo "    wrote: $out"
}

list_snapshots() {
    local dir
    if [ -n "$run_dir" ]; then
        dir="$run_dir/fields"
    else
        dir="$(pwd)/fields"
    fi
    if [ -d "$dir" ]; then
        ls "$dir"/field*.fld 2>/dev/null || true
    fi
    if [ -n "$run_dir" ] && [ -f "$run_dir/Dati.cart.out" ]; then
        printf '%s\n' "$run_dir/Dati.cart.out"
    fi
}

cmd_all() {
    echo "==> full postprocessing for run: ${run_dir:-$(pwd)}"

    # 1. statistics
    cmd_stats "$(resolve_out 'statistics.png')"

    # 2. slice each snapshot
    if [ -n "$SLICE_INDEXES" ]; then
        for i in $SLICE_INDEXES; do
            cmd_slice "$run_dir/fields/field$i.fld"
        done
    else
        while IFS= read -r f; do [ -n "$f" ] && cmd_slice "$f"; done \
            < <(list_snapshots)
    fi
}

cmd_list() {
    echo "$SEP"
    echo "==> available pipeline tools"
    "$PYTHON_BIN" -m post.python --list
}

# --- dispatch ----------------------------------------------------------------
case "$command" in
    stats)
        cmd_stats "$(resolve_out 'statistics.png')"
        ;;
    slice)
        [ $# -ge 1 ] || { echo "slice needs a <field> (e.g. fields/field1.fld)" >&2; exit 2; }
        cmd_slice "$@"
        ;;
    all)
        cmd_all
        ;;
    list)
        cmd_list
        ;;
    help|-h|--help)
        show_help
        ;;
    *)
        echo "unknown command '$command' (see $0 -h)" >&2
        exit 2
        ;;
esac
