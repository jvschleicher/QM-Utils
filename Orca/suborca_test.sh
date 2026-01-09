#!/bin/bash

# Defaults
NTASKS=1
CPUS=1
MEM=4000

# Parse options
while getopts "p:m:" opt; do
  case $opt in
    p) CPUS="$OPTARG" ;;
    m) MEM="$OPTARG" ;;
    *) echo "Usage: suborca -p <cpus> -m <mem_MB> input.inp"
       exit 1 ;;
  esac
done

shift $((OPTIND - 1))

INP="$1"

if [[ -z "$INP" || ! -f "$INP" ]]; then
  echo "Error: input file not found"
  exit 1
fi

BASENAME=$(basename "$INP" .inp)
WORKDIR=$(pwd)

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=orca_${BASENAME}
#SBATCH --output=${BASENAME}_%j.out
#SBATCH --error=${BASENAME}_%j.err
#SBATCH --ntasks=${NTASKS}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}M
#SBATCH --time=2-00:00:00
#SBATCH --signal=B:TERM@60


# module load orca

ORIG_DIR="${WORKDIR}"
SLURM_TMPDIR=/scratch/\$USER/\${SLURM_JOBID}.slurm

mkdir -p "\$SLURM_TMPDIR"
cd "\$SLURM_TMPDIR" || exit 1

# Always copy results back, no matter how we exit
cleanup() {
  echo "Copying files back to \$ORIG_DIR"
  rsync -a --ignore-missing-args *.out *.gbw *.prop *.hess *.uno *.qro *.loc *.xyz *.nbo "\$ORIG_DIR" 2>/dev/null || true
  cd "\$ORIG_DIR"
  rm -rf "\$SLURM_TMPDIR"
}
trap cleanup EXIT TERM INT

# Copy input
cp "\$ORIG_DIR/${INP}" .

# Run ORCA and put .out file always into home dir to keep track
stdbuf -oL orca ${BASENAME}.inp > "\$ORIG_DIR/${BASENAME}.out" #orca ${BASENAME}.inp > ${BASENAME}.out
joberror=\$?

exit \$joberror
EOF

