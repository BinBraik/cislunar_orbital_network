#!/bin/bash
#SBATCH --job-name=mintof-sweep
#SBATCH --account=cislunar
#SBATCH --partition=normal_q
#SBATCH --qos=owl_normal_base
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=700G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=braik@vt.edu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Unattended: build atlases (self-healing) + max-budget preview + full sweep.
# Safe to re-submit if it times out or gets killed — atlas_cache/ and
# atlas_sweep_checkpoint.mat both persist and are reused, not rebuilt.
# As of job 765432, all 13 families are already cached under this config,
# so a resubmit skips straight past the ~26h atlas-build phase.
#
# --qos=owl_normal_base: normal priority/billing. owl_normal_short is also
# available to this account (confirmed via `sacctmgr show assoc user=braik
# account=cislunar format=cluster,account,partition,qos -P`) and gets
# higher priority + higher caps for <=24h jobs, but at 2x the billing rate
# — kept in reserve, try owl_normal_base first.
#
# --time=1-00:00:00: shorter walltime requests backfill into scheduler gaps
# much faster than 4 days. Checkpoints after every completed sweep cell, so
# a resubmit (`sbatch scripts/run_sweep.slurm`) always picks up where it
# left off — expect this may take several resubmits to finish.
#
# --mem=700G: job 765432 (200G) OOM'd within minutes of starting the sweep
# proper, on cell 400/400 — the FIRST cell processed (sweep runs largest
# budget first by design) and also the LARGEST: DV_cap=0.300/Tmax=4.7124 is
# the base atlas's own bound, so atlas_derive_subset used to filter out
# almost nothing there, meaning a near-full-size second copy (Sall_sub) of
# the data existed alongside the already-resident Sall_base at that exact
# moment — the atlas build log for that run put total FRS rows across all
# 13 families at ~2.64 BILLION. run_overlap_dv_tmax_sweep.m now special-
# cases this: when a cell's budget equals the base atlas's own bounds, it
# aliases Sall_sub = Sall_base instead of calling atlas_derive_subset,
# which costs ~0 extra memory (MATLAB copy-on-write) instead of a near-
# full duplicate. That removes the specific spike that caused the OOM, so
# 700G (fits a standard 768G-RAM node, not just the 3 high-memory nodes)
# is a reasonable next try — still not verified against actual peak RSS,
# but the atlas cache is already warm so finding out is now cheap (minutes,
# not another day) if it isn't enough.
#
# N_WORKERS inside run_overlap_dv_tmax_sweep.m is set to 13 (not
# cpus-per-task-1): the subset-derivation and footprint-building parfor
# loops only ever use up to 13 workers (one per family) regardless of pool
# size, so a larger pool just means more idle MATLAB worker processes each
# carrying their own fixed memory overhead. This costs some parallelism in
# the later 78-pair loop in exchange for less baseline memory pressure.

module reset
module load MATLAB

cd "/home/braik/Agility test/FRS and BRS Method/Network - TOF Solution"

matlab -batch "try, run('scripts/run_overlap_mintof_maxbudget_preview.m'); catch ME, fprintf(2, '[slurm] Preview step failed (continuing to sweep anyway):\n%s\n', getReport(ME)); end; run('scripts/run_overlap_dv_tmax_sweep.m');"
