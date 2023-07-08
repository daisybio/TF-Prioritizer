#!/usr/bin/env python

import argparse
import os
import json

parser = argparse.ArgumentParser(description='Create directory json representation')

parser.add_argument('-d', '--directory', dest='directory', help='Directory', required=True)
parser.add_argument('-o', '--output', dest='output', help='Output file', required=True)
args = parser.parse_args()

directory = args.directory
output = args.output

def create_json_representation(recursive_directory):
    data = {}

    # Iterate through all files in the directory
    for filename in os.listdir(recursive_directory):
        filepath = os.path.join(recursive_directory, filename)

        # If it's a subdirectory, recursively call the function
        if os.path.isdir(filepath):
            data[filename] = create_json_representation(filepath)

        # Otherwise, add the file to the dictionary
        else:
            # Add the file information to the dictionary
            clean_name = os.path.splitext(filename)[0]
            clean_name = clean_name.replace(os.path.basename(recursive_directory), '')
            clean_name = clean_name.replace(directory, '')
            clean_name = clean_name.strip('_')
            data[clean_name] = "./assets/data/tf-groups/" + filepath

    # Return the data as a JSON string
    return data

# Call the function with a directory path
directory_json = create_json_representation(directory)

directory_json['name'] = os.path.basename(directory)

# Write the JSON string to a file
with open(output, 'w') as f:
    json.dump(directory_json, f, indent=4)
