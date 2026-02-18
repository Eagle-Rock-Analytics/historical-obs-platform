"""Compare station counts between S3 prefixes for each network.

Supports two modes:
  --stage qaqc   : Compare 3_qaqc_wx vs 3_qaqc_wx_v2 (zarr vs zarr)  [default]
  --stage clean   : Compare 2_clean_wx vs 3_qaqc_wx_v2 (netcdf vs zarr)

Prints a count summary table, then lists stations present in the source
but missing from 3_qaqc_wx_v2.
"""

import argparse
import boto3

BUCKET = "wecc-historical-wx"
V2_PREFIX = "3_qaqc_wx_v2"

STAGES = {
    "qaqc": {"prefix": "3_qaqc_wx", "ext": ".zarr", "label": "3_qaqc_wx"},
    "clean": {"prefix": "2_clean_wx", "ext": ".nc", "label": "2_clean_wx"},
}

NETWORKS = [
    "ASOSAWOS",
    "CAHYDRO",
    "CDEC",
    "CIMIS",
    "CNRFC",
    "CRN",
    "CW3E",
    "CWOP",
    "HADS",
    "HNXWFO",
    "HOLFUY",
    "HPWREN",
    "LOXWFO",
    "MAP",
    "MARITIME",
    "MTRWFO",
    "NCAWOS",
    "NDBC",
    "NOS-NWLON",
    "NOS-PORTS",
    "OtherISD",
    "RAWS",
    "SCAN",
    "SGXWFO",
    "SHASAVAL",
    "SNOTEL",
    "VCAPCD",
]


def list_zarr_names(s3, bucket, prefix):
    """Return a set of station IDs from zarr directories under prefix."""
    paginator = s3.get_paginator("list_objects_v2")
    station_ids = set()
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix, Delimiter="/"):
        for cp in page.get("CommonPrefixes", []):
            key = cp["Prefix"]
            if key.endswith(".zarr/"):
                basename = key.rstrip("/").rsplit("/", 1)[-1]
                if basename.endswith("_c.zarr"):
                    continue
                station_ids.add(basename.removesuffix(".zarr"))
    return station_ids


def list_nc_names(s3, bucket, prefix):
    """Return a set of station IDs from netcdf files under prefix."""
    paginator = s3.get_paginator("list_objects_v2")
    station_ids = set()
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            if key.endswith(".nc"):
                basename = key.rsplit("/", 1)[-1]
                station_ids.add(basename.removesuffix(".nc"))
    return station_ids


def list_station_ids(s3, bucket, prefix, ext):
    """Return a set of station IDs under prefix, dispatching by file extension."""
    if ext == ".zarr":
        return list_zarr_names(s3, bucket, prefix)
    elif ext == ".nc":
        return list_nc_names(s3, bucket, prefix)
    else:
        raise ValueError(f"Unsupported extension: {ext}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=list(STAGES.keys()),
        default="qaqc",
        help="Which source stage to compare against 3_qaqc_wx_v2 (default: qaqc)",
    )
    args = parser.parse_args()

    stage = STAGES[args.stage]
    src_prefix = stage["prefix"]
    src_ext = stage["ext"]
    src_label = stage["label"]

    s3 = boto3.client("s3")

    # --- Count summary table ---
    print(
        f"{'Network':<15} {src_label:>15} {'3_qaqc_wx_v2':>20} {'Diff (v2-src)':>15}"
    )
    print("-" * 70)

    total_src, total_v2 = 0, 0
    all_missing = {}  # network -> sorted list of missing station IDs

    for network in NETWORKS:
        src_ids = list_station_ids(
            s3, BUCKET, f"{src_prefix}/{network}/", src_ext
        )
        v2_ids = list_station_ids(s3, BUCKET, f"{V2_PREFIX}/{network}/", ".zarr")
        diff = len(v2_ids) - len(src_ids)
        total_src += len(src_ids)
        total_v2 += len(v2_ids)
        print(f"{network:<15} {len(src_ids):>15} {len(v2_ids):>20} {diff:>15}")

        missing = sorted(src_ids - v2_ids)
        if missing:
            all_missing[network] = missing

    print("-" * 70)
    print(
        f"{'TOTAL':<15} {total_src:>15} {total_v2:>20} {total_v2 - total_src:>15}"
    )

    # --- List stations in source but not in v2 ---
    if all_missing:
        total_missing = sum(len(v) for v in all_missing.values())
        print(f"\n\n{'='*70}")
        print(f"In {src_label} but NOT in 3_qaqc_wx_v2  ({total_missing} total)")
        print(f"{'='*70}")
        for network, missing in all_missing.items():
            print(f"\n{network} ({len(missing)} missing):")
            for name in missing:
                print(f"  {name}")

        # --- Export station IDs to a text file ---
        outfile = "missing_stations.dat"
        with open(outfile, "w") as f:
            for network, missing in all_missing.items():
                for name in missing:
                    f.write(f"{name}\n")
        print(f"\nExported {total_missing} station IDs to {outfile}")
    else:
        print(f"\nAll {src_label} stations are present in 3_qaqc_wx_v2.")


if __name__ == "__main__":
    main()
