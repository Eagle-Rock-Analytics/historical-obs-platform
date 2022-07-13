import gzip

# First, unzip .gz files.
# Function to unzip .gz files
with gzip.open("722101-99999-1973.gz", 'rb') as f:
    file_content = f.read()


# Next, find source column and filter out non-AWOS/ASOS stations.
## POS 28-28 should be 