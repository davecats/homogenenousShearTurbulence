# Next session: the production run's post-processing, then what next

Copy the block at the end as the opening message of the next session.

## Where things stand (end of 2026-10-05, session 13)

- **The production run is done** (FINDINGS.md "The production run at
  Re = 20000"): `examples/prod_re20000.in` (1536 x 1024 x 512, Re =
  20000) to S t = 100 on two A100 nodes in two segments of
  `jobs/horeka_prod.slurm` (jobs 5171460 and 5174814; 24 h and 27 h of
  queue for 11.85 h and 1.25 h of work), 132116 steps at 0.344 s/step,
  222 GB in `/hkfs/work/workspace/scratch/xt8786-hst/re20000` (20
  velocity and 20 pressure snapshots every 5 time units, the restart
  file; the workspace expires 2026-11-29, `ws_extend hst 30`).
  Statistics over S t = 30..100: S* = 7.4, -uv/q2 = 0.136,
  production/dissipation 1.005, Re_lambda = 143, dx/eta = 1.07, b_uu =
  +0.14.  Snapshots cost 19-29 s, restart files 9-18 s (twice 320 s in
  segment 2).
- **A bug found by the run and fixed:** `compute_cfl` never reset
  `cfl`, so the CFL-chosen step could only shrink (the run used dt =
  0.00074 where 0.00091 was allowed: 20% extra steps, always safe).
  Fixed-step decks were unaffected; safety net green; `~/hst` on HoreKA
  rebuilt with the fix.  Every earlier CFL-driven timing (the small
  runs, WP5's long run) carried the same handicap.
- **The 1e-4 between CPL's and our pressure** is the fourth-order
  y-truncation of the compact scheme (17.7x per halving of dy, in the
  highest y-modes; D1 and D2 are fourth order, D4 sixth; README
  corrected).
- **The H100 is gone for good**; HoreKA is being migrated.  Stale
  copies on HoreKA (`~/hst-y`, `~/hst-exp`, `~/hst-exp2`, `~/hst-base`,
  `~/hst-exp4`, `~/hst-exp5`) can go (`~/hst-y` holds the job outputs
  `hst-2node-517016*.out` that FINDINGS.md "Scaling" cites; copy them to
  `~/hst` first).
- **The post-processing is now Fortran** (`src/postpro/postpro.f90`, `make
  postpro`, `postpro.in`; FINDINGS.md "Post-processing in Fortran"):
  means, stresses, spectra and the six Reynolds-stress budgets per plane,
  text output, checked against the CPL chain to 5e-11 on one snapshot.
  Job 5183609 (`accelerated`, one A100 node, 40 min) runs it on the
  production run's snapshots 6..20 into `<run dir>/statistics_f90/`
  (`hst-<jobid>.out` in the run directory); the `dev_accelerated` queue
  was full of the user's jobs.
- **Pending at the end of the session:** job 5181558
  (`jobs/horeka_postprocess.slurm`, `cpuonly`, 4 ranks, `--mem=230gb`,
  6 h) runs the CPL chain on snapshots 6..20 of the production run
  (`NFMIN=6`, S t = 30..100); its output is `~/hst/hst-post-5181558.out`
  and `<run dir>/statistics/`.  It had been pending for over an hour
  when the session ended.

## The task

1. **Read the post-processing**: job 5183609 (Fortran, `statistics_f90/`:
   `stresses.dat` plane averages must give q2 about 0.083 and -uv/q2
   about 0.136, the box budgets printed at the end of its log should
   balance to a few per cent over 15 snapshots) and, as the cross-check,
   the CPL job 5181558: if it ran, `rms.dat`
   (plane averages must give q2 about 0.083 and, in its `uw` column =
   our u v, -uv/q2 about 0.136 over the 15 snapshots), `mean.dat`,
   `spectra.bin`, `uiuj.bin`, `mke.bin` (the budget: production,
   dissipation, pressure-strain per plane; its reader is in the CPL
   code's conventions, `statistics/README.md` of `hst-main` says
   nothing).  If it ran out of memory, `RANKS=2` (about 50 GB per rank
   at this size); if it timed out, split by `NFMIN`/`NFMAX`.  Put the
   numbers next to the `Runtimedata` ones in FINDINGS.md.
2. **Decide what the next run is.**  Candidates: (a) the same deck
   again with the CFL fix (20% fewer steps, one 12 h segment should
   reach S t = 100); (b) a higher Re on the same grid is not resolved
   (dx/eta = 1.07 already), so a higher Re_lambda needs a finer grid
   with the same shape: nx = 3 nz with ny = 2 (nx + 1)/3... in practice
   `nx = 767, ny = 1536, nz = 255` (dx = dy = dz = 0.0013, 3.4x the
   points, about 44 GB per GPU on 8 ranks: four nodes, or two with nx =
   nz = 511 and ny = 1024 as the 1024^3 benchmark but z over-resolved);
   eta ~ Re^(-3/4) says Re = 35000 for dx/eta = 1.1 on that grid; (c) the
   S2 or Stokes-layer variants of the deck (`&physics`), whose
   post-processing the CPL chain also knows.
3. **Later items unchanged:** the reduced-system gather for `npy > 2`
   (FINDINGS.md "Scaling": 18% of the four-node 1024^3 step), the I/O
   (one write of 8.6 GB took 320 s twice: MPI-IO hints), the cleanup of
   the HoreKA copies, the repeated `Runtimedata` line at a restart.

## Safety net (the third argument is npy)

```bash
tests/run_tests.sh build-cpu 2 && tests/run_tests.sh build-cpu 4 2 && tests/run_tests.sh build-cpu 4 4   # 12 runs each
tests/run_tests.sh build-gpu 2 && tests/run_tests.sh build-gpu 2 2
tests/regression.sh build-cpu 2 && tests/regression.sh build-cpu 4 2        # 50 steps of three decks at 1e-10
tests/regression.sh build-gpu 1 && tests/regression.sh build-gpu 2 && tests/regression.sh build-gpu 2 2
tests/crossbranch.sh build-cpu build-cpu 4 1 2 && tests/crossbranch.sh build-cpu build-cpu 4 2 1   # restart between npy values
ssh istmcetus 'cd ~/Codes/hst/homogenenousShearTurbulence && source env/istm.sh && make GPU=1 NCCL=1 BUILD=build-nccl -j8 && make GPU=1 NCCL=1 BUILD=build-nccl test -j8 && tests/run_tests.sh build-nccl 2 2 && tests/regression.sh build-nccl 2 2 && tests/run_tests.sh build-nccl 2 && tests/regression.sh build-nccl 2 && tests/crossbranch.sh build-nccl build-nccl 2 1 2'
```

The CPU build must be made in a shell that has *not* sourced
`env/istm.sh`, the GPU tests in one that *has*; `make test` rebuilds the
test programs; a flag change needs the objects removed.  Parallel Bash
calls share one working directory: absolute paths or `make -C`.  A
bit-identity check of a change that must not alter the numbers: the
`small` deck with the GPU build before and after (50 steps), `cmp` the
two `Dati.cart.out`.  nvfortran 25.9 rejects the names `kind` and `x` in
OpenMP clauses of any file that can see `hst_mpi`; an absent optional
dummy must not appear in a target region.  HoreKA: `rsync -a --delete
--exclude 'build-*' --exclude .git --exclude 'hst-*.out' ./ horeka:hst/`,
then `source env/horeka.sh gpu; make GPU=1 GPU_ARCH=cc80 NCCL=1`; never
rebuild a root while a job is queued on it (the production chain reads
`~/hst/build-gpu/hst` at every segment: rebuild `~/hst` only between
segments, or point the chain at a copy with `HST_ROOT`).  Scripted ssh to
HoreKA: `ssh -o BatchMode=yes -o ProxyCommand=false horeka ...`, one call
at a time; `ssh -O check horeka` says whether the master is alive; when
it is not, the user types `! ssh horeka true`.

## Opening message for the next session

```
Repository ~/Codes/hst/homogenenousShearTurbulence (one branch, main; on
HoreKA ~/hst = main), a GPU/CPU DNS for homogeneous shear turbulence with
x-z pencils times npy y slabs; read README.md, then NEXT_SESSION.md, then
FINDINGS.md from "The production run at Re = 20000" to the end.  Task:
NEXT_SESSION.md item 1 (the CPL post-processing of the production run),
then item 2 (decide the next run with me) unless I say otherwise.  Safety net
green after any code change (CPU, GPU, NCCL on istmcetus, the npy restart
round trip).  Do not modify ~/Codes/hst/channel or ~/Codes/hst/hst-main.
Commit each step; push at the end.  First thing: `! ssh horeka true` in
the prompt if `ssh -O check horeka` says the master is gone.
```
