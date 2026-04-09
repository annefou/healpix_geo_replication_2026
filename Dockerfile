# HEALPix-geo Benchmark Replication Environment
#
# Replication of: Law & Ardo (2024)
# DOI: 10.1080/20964471.2024.2429847
# Extension: HEALPix indexing via healpix-geo (sphere + WGS84 ellipsoid)
#
# Build:  docker build -t healpix-geo-replication .
# Run:    docker run -v $(pwd)/results:/app/results healpix-geo-replication
#
# Environment variables (passed via docker run -e):
#   VECTOR_LAYERS     Comma-separated list (e.g., "5,10,20,50,100")
#   RASTER_LAYERS     Comma-separated list (e.g., "10,50,100,500")
#   HEALPIX_DEPTH     HEALPix depth (default: 9)
#   RANDOM_SEED       Random seed for reproducibility (default: 42)
#
# Author: Anne Fouilloux
# Date: 2026-04-09

FROM python:3.11-slim-bookworm

# Set labels
LABEL maintainer="Anne Fouilloux <annefou@lifewatch.eu>"
LABEL description="Reproducible HEALPix-geo benchmark replication (sphere + WGS84)"
LABEL paper.doi="10.1080/20964471.2024.2429847"

# Create app directory
WORKDIR /app

# Install system dependencies for geopandas/shapely
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgeos-dev \
    libproj-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY run_healpix_geo_replication.py /app/
COPY run_comparison.py /app/
COPY config.env /app/
COPY entrypoint.sh /app/
COPY reference_h3/ /app/reference_h3/
COPY README.md /app/

# Create directories for results
RUN mkdir -p /app/results /app/results_healpix_geo /app/results_comparison

# Set environment variables for reproducibility
ENV PYTHONHASHSEED=42
ENV PYTHONUNBUFFERED=1

# Default configuration (can be overridden via docker run -e)
ENV VECTOR_LAYERS="5,10,20,50,100"
ENV RASTER_LAYERS="10,50,100,500"
ENV HEALPIX_DEPTH="9"
ENV RANDOM_SEED="42"

# Make entrypoint executable
RUN chmod +x /app/entrypoint.sh
ENTRYPOINT ["/app/entrypoint.sh"]

# Default: run healpix-geo benchmarks + comparison
CMD ["sh", "-c", "\
    python run_healpix_geo_replication.py --all --output results && \
    python run_comparison.py --healpix-geo results --output results \
"]
