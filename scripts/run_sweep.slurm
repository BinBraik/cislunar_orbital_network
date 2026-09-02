#!/bin/bash
#SBATCH --job-name=mintof-sweep
#SBATCH --account=cislunar
#SBATCH --partition=normal_q
#SBATCH --qos=owl_normal_short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=3000G
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
# --qos=owl_normal_short: confirmed available to this account on owl
# (`sacctmgr show assoc user=braik account=cislunar format=cluster,account,
# partition,qos -P`) — higher priority + higher resource caps than
# owl_normal_base, specifically for <=24h jobs, matching --time below.
#
# --time=1-00:00:00: shorter walltime requests backfill into scheduler gaps
# much faster than 4 days. Checkpoints after every completed sweep cell, so
# a resubmit (`sbatch scripts/run_sweep.slurm`) always picks up where it
# left off — expect this may take several resubmits to finish.
#
# --mem=3000G: job 765432 (200G) OOM'd within minutes of starting the sweep
# proper, on cell 400/400 — the FIRST cell processed (sweep runs largest
# budget first by design) and also the LARGEST: DV_cap=0.300/Tmax=4.7124 is
# the base atlas's own bound, so atlas_derive_subset filters out almost
# nothing there, meaning a near-full-size second copy (Sall_sub) of the
# data exists alongside the already-resident Sall_base at that exact
# moment. The atlas build log for that run put total FRS rows across all
# 13 families at ~2.64 BILLION — far larger than assumed when 1200G was
# picked (borrowed from a different, unrelated job). 3000G is a safer
# margin given that scale, not a value verified against actual peak RSS.
#
# N_WORKERS inside run_overlap_dv_tmax_sweep.m is set to 13 (not
# cpus-per-task-1): the subset-derivation and footprint-building parfor
# loops only ever use up to 13 workers (one per family) regardless of pool
# size, so a larger pool just means more idle MATLAB worker processes each
# carrying their own fixed memory overhead during exactly the phase that
# OOM'd. This costs some parallelism in the later 78-pair loop in exchange
# for less baseline memory pressure at the riskiest point.

module reset
module load MATLAB

cd "/home/braik/Agility test/FRS and BRS Method/Network - TOF Solution"

matlab -batch "try, run('scripts/run_overlap_mintof_maxbudget_preview.m'); catch ME, fprintf(2, '[slurm] Preview step failed (continuing to sweep anyway):\n%s\n', getReport(ME)); end; run('scripts/run_overlap_dv_tmax_sweep.m');"
