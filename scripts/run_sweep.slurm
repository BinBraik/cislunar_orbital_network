#!/bin/bash
#SBATCH --job-name=mintof-sweep
#SBATCH --partition=owl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=mintof-sweep_%j.out
#SBATCH --error=mintof-sweep_%j.err
# Unattended: build atlases (self-healing) + max-budget preview + full sweep.
# Safe to re-submit if it times out or gets killed — atlas_cache/ and
# atlas_sweep_checkpoint.mat both persist and are reused, not rebuilt.
#
# BEFORE SUBMITTING, CHECK / ADJUST:
#   --partition   : confirm "owl" is correct for your allocation
#   --cpus-per-task, --mem, --time : verified against `sinfo -p owl` limits,
#       not against actual atlas file sizes (unknown from here) — if this
#       OOMs or times out, bump --mem / --time and resubmit; it resumes.
#   module load matlab : exact module name/version — check `module avail matlab`
#       on Owl and edit the line below if it differs.
#   N_WORKERS inside run_overlap_dv_tmax_sweep.m is set to 7 to match
#       --cpus-per-task=8 (7 parfor workers + 1 main thread). If you change
#       --cpus-per-task, change N_WORKERS to match (leave 1 core spare).

set -euo pipefail

# Repo root — this script lives in <repo>/scripts/, so go one level up.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

module load matlab

echo "[slurm] Job $SLURM_JOB_ID starting in: $REPO_DIR"
echo "[slurm] Nodes: $SLURM_JOB_NODELIST   CPUs/task: $SLURM_CPUS_PER_TASK   Mem: ${SLURM_MEM_PER_NODE:-unset}"

matlab -nodisplay -nosplash -nodesktop -r "\
try, \
    run('scripts/run_overlap_mintof_maxbudget_preview.m'); \
catch ME, \
    fprintf(2, '[slurm] Preview step failed (continuing to sweep anyway):\n%s\n', getReport(ME)); \
end; \
run('scripts/run_overlap_dv_tmax_sweep.m'); \
exit;"

echo "[slurm] Job $SLURM_JOB_ID finished."
