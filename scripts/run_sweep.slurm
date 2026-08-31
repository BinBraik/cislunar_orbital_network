#!/bin/bash
#SBATCH --job-name=mintof-sweep
#SBATCH --account=cislunar
#SBATCH --partition=normal_q
#SBATCH --qos=owl_normal_base
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=200G
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=braik@vt.edu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Unattended: build atlases (self-healing) + max-budget preview + full sweep.
# Safe to re-submit if it times out or gets killed — atlas_cache/ and
# atlas_sweep_checkpoint.mat both persist and are reused, not rebuilt.
# N_WORKERS inside run_overlap_dv_tmax_sweep.m is set to 39 to match
# --cpus-per-task=40 (39 parfor workers + 1 main thread). If you change
# --cpus-per-task, change N_WORKERS to match (leave 1 core spare).

module reset
module load MATLAB

cd "/home/braik/Agility test/FRS and BRS Method/Network - TOF Solution"

matlab -batch "try, run('scripts/run_overlap_mintof_maxbudget_preview.m'); catch ME, fprintf(2, '[slurm] Preview step failed (continuing to sweep anyway):\n%s\n', getReport(ME)); end; run('scripts/run_overlap_dv_tmax_sweep.m');"
