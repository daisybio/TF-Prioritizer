#!/bin/sh

touch {STATSFILE}
docker stats --no-stream | sed -n 1p | sed 's/ \{2,\}/\t/g' | sed 's/^/TIME\t/' | cat >> {STATSFILE}
while true; do
  sleep 5;
  docker stats --no-stream | sed -n 2p | sed 's/ \{2,\}/\t/g' | sed "s/^/$(date)\t/" | cat >> {STATSFILE}
done