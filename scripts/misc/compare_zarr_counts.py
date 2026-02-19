"""Compare station counts between two S3 prefixes.

Usage:
    python compare_zarr_counts.py <prefix1> <prefix2> <type1> <type2>

Arguments:
    prefix1   First S3 prefix (e.g. 3_qaqc_wx)
    prefix2   Second S3 prefix (e.g. 3_qaqc_wx_v2)
    type1     Data type for prefix1: zarr or nc
    type2     Data type for prefix2: zarr or nc

Examples:
    # Compare two qaqc zarr prefixes
    python compare_zarr_counts.py 3_qaqc_wx 3_qaqc_wx_v2 zarr zarr

    # Compare clean netcdf against qaqc zarr
    python compare_zarr_counts.py 2_clean_wx 3_qaqc_wx_v2 nc zarr

    # Same as above, also export missing stations to missing_stations.dat
    python compare_zarr_counts.py 2_clean_wx 3_qaqc_wx_v2 nc zarr --export

    # Compare qaqc zarr against merge zarr, export missing stations to .dat
    python compare_zarr_counts.py 3_qaqc_wx_v2 4_merge_wx_v2 zarr zarr --export

    # Use a different S3 bucket
    python compare_zarr_counts.py 3_qaqc_wx 3_qaqc_wx_v2 zarr zarr --bucket my-other-bucket
"""

import argparse
import boto3

BUCKET = "wecc-historical-wx"

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


def list_station_ids(s3, bucket, prefix, dtype):
    if dtype == "zarr":
        return list_zarr_names(s3, bucket, prefix)
    elif dtype == "nc":
        return list_nc_names(s3, bucket, prefix)
    else:
        raise ValueError(f"Unsupported type: {dtype}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("prefix1", help="First S3 prefix (e.g. 3_qaqc_wx)")
    parser.add_argument("prefix2", help="Second S3 prefix (e.g. 3_qaqc_wx_v2)")
    parser.add_argument("type1", choices=["zarr", "nc"], help="Data type for prefix1")
    parser.add_argument("type2", choices=["zarr", "nc"], help="Data type for prefix2")
    parser.add_argument(
        "--bucket", default=BUCKET, help=f"S3 bucket name (default: {BUCKET})"
    )
    parser.add_argument(
        "--export",
        action="store_true",
        default=False,
        help="Export missing station IDs to missing_stations.dat",
    )
    args = parser.parse_args()

    s3 = boto3.client("s3")

    col1 = args.prefix1[-20:]  # truncate for display
    col2 = args.prefix2[-20:]
    print(f"{'Network':<15} {col1:>20} {col2:>20} {'Diff (2-1)':>12}")
    print("-" * 70)

    total1, total2 = 0, 0
    all_missing = {}

    for network in NETWORKS:
        ids1 = list_station_ids(
            s3, args.bucket, f"{args.prefix1}/{network}/", args.type1
        )
        ids2 = list_station_ids(
            s3, args.bucket, f"{args.prefix2}/{network}/", args.type2
        )
        diff = len(ids2) - len(ids1)
        total1 += len(ids1)
        total2 += len(ids2)
        print(f"{network:<15} {len(ids1):>20} {len(ids2):>20} {diff:>12}")

        missing = sorted(ids1 - ids2)
        if missing:
            all_missing[network] = missing

    print("-" * 70)
    print(f"{'TOTAL':<15} {total1:>20} {total2:>20} {total2 - total1:>12}")

    if all_missing:
        total_missing = sum(len(v) for v in all_missing.values())
        print(f"\n\n{'='*70}")
        print(f"In {args.prefix1} but NOT in {args.prefix2}  ({total_missing} total)")
        print(f"{'='*70}")
        for network, missing in all_missing.items():
            print(f"\n{network} ({len(missing)} missing):")
            for name in missing:
                print(f"  {name}")

        if args.export:
            outfile = "missing_stations.dat"
            with open(outfile, "w") as f:
                for network, missing in all_missing.items():
                    for name in missing:
                        f.write(f"{name}\n")
            print(f"\nExported {total_missing} station IDs to {outfile}")
    else:
        print(f"\nAll {args.prefix1} stations are present in {args.prefix2}.")


if __name__ == "__main__":
    main()
