#!/usr/bin/env python3
"""
Download CAMI marine mock community dataset.
If URL is provided, download and extract; otherwise create a placeholder.
"""

import os
import sys
import argparse
import urllib.request
import tarfile


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    p.add_argument("--tmpdir", required=True)
    p.add_argument("--url", default="")
    args = p.parse_args()

    if args.url:
        try:
            archive = os.path.join(args.tmpdir, "cami_marine.tar.gz")
            print(f"Downloading CAMI from {args.url} ...")
            urllib.request.urlretrieve(args.url, archive)
            print("Extracting ...")
            with tarfile.open(archive, "r:gz") as tar:
                tar.extractall(path=args.tmpdir)
            # Move contents to outdir
            for item in os.listdir(args.tmpdir):
                if item != "cami_marine.tar.gz":
                    src = os.path.join(args.tmpdir, item)
                    dst = os.path.join(args.outdir, item)
                    if os.path.exists(dst):
                        # merge or replace
                        import shutil

                        shutil.rmtree(dst) if os.path.isdir(dst) else os.remove(dst)
                    shutil.move(src, dst)
            print("CAMI data ready.")
            return
        except Exception as e:
            print(f"Failed to download/extract CAMI: {e}")
            # Fall through to create placeholder

    # Placeholder
    readme = os.path.join(args.outdir, "README.txt")
    with open(readme, "w") as f:
        f.write("CAMI Marine mock community placeholder.\n")
        f.write(
            "To use real data, download from https://data.cami-challenge.org/ and place here.\n"
        )
    print("CAMI placeholder created.")


if __name__ == "__main__":
    main()
