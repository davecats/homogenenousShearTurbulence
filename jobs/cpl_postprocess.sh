#!/bin/bash
# Post-process a run directory with the CPL tools of hst-main (postpro.cpl:
# Reynolds stresses and spectra per y plane, energy budgets, optionally the
# pressure reconstruction), which read our files unchanged.
#   jobs/cpl_postprocess.sh <run dir> [nranks] [nfmin nfmax dn]
# Writes scddns.in (the deck in CPL names: their ny = our nz, their nz = our
# ny + 1) next to hst.in and CPL's postpro.in in <run dir>/cpl-postpro (our own
# postpro.in, the namelist of src/postpro, keeps its name), builds postpro there
# from a copy of HST_MAIN (~/Codes/hst/hst-main; ~/hst-main on HoreKA) with
# mpicpl (needs cpl and mpicc on PATH) and runs it with mpirun; the results go
# to <run dir>/statistics/ (mean.dat, rms.dat, spectra.bin, uiuj.bin, mke.bin).
# The pressure snapshots p_fields/ written by the run are used as they are;
# PRESSURE=cpl recomputes them with prepare_pressure.cpl instead (overwriting
# ours: copy them aside first).  The copy of hst-main gets one patch: its
# postprocess/ and pressure_reconstruction/ call penta_smw_solve with a fifth
# argument that its linsolver_smw.cpl no longer takes.
set -eu
run=$(cd "$1" && pwd); np=${2:-1}
HST_MAIN=${HST_MAIN:-$HOME/Codes/hst/hst-main}
[ -f "$run/hst.in" ] || { echo "no hst.in in $run"; exit 1; }

# the deck's values (namelist: name = value, comments after !), with hst_input's defaults
val() { v=$(sed 's/!.*//' "$run/hst.in" | grep -oE "\b$1\b *= *[^,/ ]+" | tail -1 | sed 's/.*= *//'); echo "${v:-$2}"; }
nx=$(val nx 63); ny=$(val ny 128); nz=$(val nz 63)
alfa0=$(val alfa0 2.0943951023931953); beta0=$(val beta0 6.283185307179586); ystretch=$(val ystretch 0.0)
re=$(val re 1000.0); S=$(val S 1.0); s2a=$(val s2_amplitude 0.0); s2p=$(val s2_period 0.0); s2s=$(val s2_start 0.0)
deltat=$(val deltat 0.0); cflmax=$(val cflmax 1.0); t_max=$(val t_max 100.0)
dt_stat=$(val dt_stat 0.01); dt_field=$(val dt_field 10.0); dt_save=$(val dt_save 10.0)
htcoeff=-1; awk "BEGIN {exit !($ystretch > 0)}" && htcoeff=$ystretch
cat > "$run/scddns.in" <<END
nx=$nx   ny=$nz   nz=$((ny + 1))
alpha0=$alfa0 beta0=$beta0
htcoeff=$htcoeff
Re=$re Pr=0.71
deltat=$deltat   cflmax=$cflmax	t_max=$t_max
u_conv=0 v_conv=0
dt_field=$dt_field    dt_save=$dt_save dt_stat=$dt_stat   dt_chk_cfl=0.01	time_from_restart=NO
u0=0 un=0 v0=0 vn=0 phi0=-1 phin=1
S=$S A=$s2a T=$s2p delta_SL=0.02 t0_SL=$s2s
meanpx=0	meanpy=0
Vfield=
Sfield=
END
nf=$(ls "$run"/fields/field*.fld 2>/dev/null | wc -l)
[ "$nf" -gt 0 ] || { echo "no fields/field*.fld in $run"; exit 1; }
b="$run/cpl-postpro"; mkdir -p "$b"
printf 'nfmin=%s\nnfmax=%s\ndn=%s\npath_name=%s/\n' "${3:-1}" "${4:-$nf}" "${5:-1}" "$run" > "$b/postpro.in"
echo "$run: $nf fields, CPL deck nx=$nx ny=$nz nz=$((ny + 1)), $(tr '\n' ' ' < "$b/postpro.in")"

# build postpro from a copy of hst-main (the pressure flag off unless PRESSURE=cpl)
if [ ! -x "$b/postpro" ] || [ "${PRESSURE:-}" != "$(cat "$b/.pressure" 2>/dev/null)" ]; then
  rm -rf "$b"/*.cpl "$b"/postprocess "$b"/pressure_reconstruction "$b"/.cpl "$b"/postpro
  cp -r "$HST_MAIN"/*.cpl "$HST_MAIN"/postprocess "$HST_MAIN"/pressure_reconstruction "$b"
  chmod -R u+w "$b"
  sed -i 's/, check_linsolve)/)/' "$b"/postprocess/convenience.cpl "$b"/pressure_reconstruction/poisson.cpl
  [ "${PRESSURE:-}" = cpl ] || sed -i 's/^#define pressure_fields/! #define pressure_fields/' "$b/flags.cpl"
  echo "${PRESSURE:-}" > "$b/.pressure"
  # mpicpl compiles with the flags of `mpicc -show` but links with -lmpi alone: give the linker the -L paths (HoreKA's modules)
  ( cd "$b" && LOADLIBES="$(mpicc -show | grep -oE -- '-L[^ ]+' | tr '\n' ' ')${LOADLIBES:-}" mpicpl postpro.cpl 2>&1 \
    | sed 's/\x1b\[[0-9;?]*[a-zA-Z]//g' | grep -E "ERROR|error" ) || true
  [ -x "$b/postpro" ] || { echo "postpro did not build (see $b)"; exit 1; }
fi
mkdir -p "$run/statistics"
( cd "$b" && mpirun -np "$np" ./postpro )
ls -la "$run/statistics"
