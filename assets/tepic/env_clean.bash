#!/bin/bash


# Remove version annotations from environment.yml
sed -i 's/=.*$/=/' environment.yml

echo "Version annotations removed from environment.yml"
