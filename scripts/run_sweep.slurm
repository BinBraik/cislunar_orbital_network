#!/bin/bash
#SBATCH --job-name=mintof-sweep
#SBATCH --account=cislunar
#SBATCH --partition=normal_q
#SBATCH --qos=owl_normal_base
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=1200G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=braik@vt.edu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Unattended: build atlases (self-healing) + max-budget preview + full sweep.
# Safe to re-submit if it times out or gets killed — atlas_cache/ and
# atlas_sweep_checkpoint.mat both persist and are reused, not rebuilt.
#
# --time=1-00:00:00 (shortened from 4 days): shorter walltime requests
# backfill into scheduler gaps much faster. The 200G attempt (job 765432)
# already took 1d2h before OOMing, so this WILL very likely hit the 1-day
# limit before finishing the full sweep — that's fine, it checkpoints
# after every completed cell, so just resubmit (`sbatch scripts/run_sweep.slurm`)
# and it picks up where it left off. Expect several resubmits to finish.
# N_WORKERS inside run_overlap_dv_tmax_sweep.m is set to 39 to match
# --cpus-per-task=40 (39 parfor workers + 1 main thread). If you change
# --cpus-per-task, change N_WORKERS to match (leave 1 core spare).
#
# --mem=1200G (raised back from a 200G attempt that OOM'd after ~1 day —
# job 765432, state out_of_memory). Sall_base holds all 13 full base
# atlases in memory for the WHOLE run, and at this extended budget
# (0.3/3*pi/2) each atlas logs substantially more hit rows than the
# paper's nominal 0.2/pi config, so memory scales with atlas size/worker
# count, not core count. 1200G matches this account's own previously
# working large-memory job (run_ta_stable_stable_pairs_extended.slurm).

module reset
module load MATLAB

cd "/home/braik/Agility test/FRS and BRS Method/Network - TOF Solution"

matlab -batch "try, run('scripts/run_overlap_mintof_maxbudget_preview.m'); catch ME, fprintf(2, '[slurm] Preview step failed (continuing to sweep anyway):\n%s\n', getReport(ME)); end; run('scripts/run_overlap_dv_tmax_sweep.m');"
