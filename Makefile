build:
	docker build -t af2_msa:alpha -f docker/Dockerfile . \
	--build-arg USER_ID=$(id -u) \
	--build-arg GROUP_ID=$(id -g)

run:
	docker run -it --rm --name=af2_msa \
	-v $(PWD):/app/alphafold \
	-v $(PWD)/test:/app/data \
	-v $(PWD)/params:/app/alphafold/alphafold/data/params/ \
	af2_msa:alpha \
	python3 /app/alphafold/run_alphafold.py \
	--precomputed_msa /app/data/PF02518_seed.fasta \
	--output_dir /app/data/output \
	--data_dir /app/data/data
