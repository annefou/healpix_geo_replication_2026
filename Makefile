# HEALPix-geo Benchmark Replication - Makefile
#
# Usage:
#   make help          Show available commands
#   make docker-build  Build Docker image
#   make docker-run    Run replication in Docker
#   make local-setup   Setup local Python environment
#   make run           Run replication locally
#   make clean         Clean generated files

.PHONY: help docker-build docker-run docker-shell local-setup run run-benchmark run-comparison clean clean-results

# Load configuration from config.env
include config.env
export

# Default target
help:
	@echo "HEALPix-geo Benchmark Replication Study"
	@echo "========================================"
	@echo ""
	@echo "Docker commands:"
	@echo "  make docker-build   Build Docker image"
	@echo "  make docker-run     Run benchmarks + comparison in Docker"
	@echo "  make docker-shell   Open shell in Docker container"
	@echo ""
	@echo "Local commands:"
	@echo "  make local-setup    Create Python venv and install deps"
	@echo "  make run            Run benchmarks + comparison locally"
	@echo "  make run-benchmark  Run healpix-geo benchmarks only"
	@echo "  make run-comparison Run comparison against H3 reference only"
	@echo ""
	@echo "Utility commands:"
	@echo "  make clean          Remove generated files"
	@echo "  make show-config    Show current configuration"
	@echo ""

# Show configuration
show-config:
	@echo "Current configuration (from config.env):"
	@echo "  HEALPIX_DEPTH:  $(HEALPIX_DEPTH)"
	@echo "  VECTOR_LAYERS:  $(VECTOR_LAYERS)"
	@echo "  RASTER_LAYERS:  $(RASTER_LAYERS)"
	@echo "  RANDOM_SEED:    $(RANDOM_SEED)"

# Docker commands
DOCKER_IMAGE = healpix-geo-replication
DOCKER_TAG = latest

docker-build:
	docker build -t $(DOCKER_IMAGE):$(DOCKER_TAG) .

docker-run: docker-build
	docker run --rm \
		-v $(PWD)/results:/app/results \
		-e RANDOM_SEED=$(RANDOM_SEED) \
		$(DOCKER_IMAGE):$(DOCKER_TAG)

docker-shell: docker-build
	docker run -it --rm \
		-v $(PWD)/results:/app/results \
		$(DOCKER_IMAGE):$(DOCKER_TAG) \
		bash

# Local Python environment
VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip

$(VENV)/bin/activate: requirements.txt
	python3 -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt
	touch $(VENV)/bin/activate

local-setup: $(VENV)/bin/activate
	@echo "Virtual environment ready. Activate with: source $(VENV)/bin/activate"

# Run commands (local)
run: $(VENV)/bin/activate
	$(PYTHON) run_healpix_geo_replication.py --all --output results
	$(PYTHON) run_comparison.py --healpix-geo results --output results

run-benchmark: $(VENV)/bin/activate
	$(PYTHON) run_healpix_geo_replication.py --all --output results

run-comparison: $(VENV)/bin/activate
	$(PYTHON) run_comparison.py --healpix-geo results --output results

# Cleanup
clean:
	rm -rf results/
	rm -rf results_healpix_geo/
	rm -rf results_comparison/
	rm -rf $(VENV)
	rm -rf __pycache__
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete

clean-results:
	rm -rf results/
	rm -rf results_healpix_geo/
	rm -rf results_comparison/
