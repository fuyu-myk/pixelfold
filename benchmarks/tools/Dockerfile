# Pinned reference tools for pixelfold golden generation: DSSP 4 + FreeSASA.
#
# Build (from the repo root):
#   docker build -t pixelfold-golden .
# Run:
#   docker run --rm -v "$PWD/benchmarks:/work/benchmarks" pixelfold-golden \
#       --manifest benchmarks/manifest.toml --golden-dir benchmarks/golden
#
# Pin versions and platforms here so a regenerated golden set is reproducible. If a pin fails
# to resolve, bump it and record the change alongside the regenerated files.
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8

# mkdssp (DSSP 4), Python, Biopython (DSSP parsing) and FreeSASA (SASA).
RUN micromamba install -y -n base -c conda-forge \
        python=3.11 \
        dssp=4.5.8 \
        biopython=1.84 \
        freesasa-python=2.2.1 \
    && micromamba clean --all --yes

COPY generate_golden.py /usr/local/bin/generate_golden.py

WORKDIR /work
# Route through micromamba's _entrypoint.sh so the base env is activated.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/usr/local/bin/generate_golden.py"]
