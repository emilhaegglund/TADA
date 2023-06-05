# Script to symlink files from the NCBI and the Prokka directories
import os
import sys

print(snakemake.output[0])
out_dir = os.path.abspath(snakemake.output[0])

# Make sure that the output dir exists
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
for f in snakemake.input:
    fname = os.path.split(f)[1]
    src_abs_path = os.path.abspath(f)
    target = os.path.join(out_dir,fname)
    os.symlink(src=src_abs_path, dst=target, target_is_directory=False)
    print(f"Symlinking {src_abs_path} to {out_dir}")

