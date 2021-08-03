Singularity is a common container tool on HPC environments which often don't allow docker.
This script should work with any Singularity installation that is newer than version 3.2.

# Example build 
singularity build alphafold.simg alphafold.singularity

# Example run
./alphafold.simg python /app/alphafold/run_alphafold.py
