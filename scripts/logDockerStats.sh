#!/bin/sh

touch {STATSFILE}
echo "TIME\tCONTAINER_ID\tCONTAINER_NAME\tCPU_PERCENTAGE\tMEMORY_USAGE\tMEMORY_PERCENTAGE\tPIDs" >> {STATSFILE}
while true; do
  sleep 5;
  docker stats --no-stream --format "table {{.ID}}\t{{.Name}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}\t{{.PIDs}}" \
  {CONTAINER_ID} | sed -n 2p | sed "s/^/$(date)\t/" | cat >> {STATSFILE}
done