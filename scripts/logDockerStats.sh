#!/bin/sh

touch {STATSFILE}
echo "TIME\tCONTAINER_ID\tCONTAINER_NAME\tCPU_PERCENTAGE\tMEMORY_USAGE\tMEMORY_PERCENTAGE\tPIDs" >> {STATSFILE}
while true; do
  sleep 5;
  docker stats --no-stream {CONTAINER_ID} | sed -n 2p | sed "s/^/$(date)\t/" | cat >> {STATSFILE}
done