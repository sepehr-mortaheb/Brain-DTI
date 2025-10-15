
BIDS_DIR="./data/"
OUT_DIR="./qc_output/"
WORK_DIR="./scratch/mriqc_work/"
TF_DIR="$HOME/.cache/templateflow/"

docker run --rm -it \
  -v "$BIDS_DIR":/data:ro \
  -v "$OUT_DIR":/out \
  -v "$WORK_DIR":/work \
  -v "$TF_DIR":/opt/templateflow \
  -e TEMPLATEFLOW_HOME=/opt/templateflow \
  -u $(id -u):$(id -g) \
  nipreps/mriqc:24.0.2 \
  /data /out participant \
  -m T1w bold dwi \
  -w /work \
  --participant-label cosmonaut03 \
  --no-sub

