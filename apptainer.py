import argparse
import json
import os

parser = argparse.ArgumentParser()

parser.add_argument('-c', type=str, help='Path to config file')
parser.add_argument('-t', type=str, help='Thread count')
parser.add_argument('-o', type=str, help='Output directory')
parser.add_argument('-m', type=int, default=10, help='Memory limit in GB, default: 10')
parser.add_argument('-i', type=str,
                    help='[Optional] Docker image, only for development purposes')
parser.add_argument('-s', type=str, default=None, help='Slurm partition')

args = parser.parse_args()


def process_configs(configs: dict) -> dict:
    mounts = {}
    for key, value in configs.items():
        if isinstance(value, dict):
            mounts.update(process_configs(value))
        else:
            if isinstance(value, str) and os.path.exists(value):
                extension = os.path.splitext(value)[1]
                inside_path = f"/srv/input/{key}{extension}"
                mounts[inside_path] = os.path.abspath(value)

                configs[key] = inside_path
    return mounts


image = args.i if args.i else "nicotru/tf-prioritizer:master"
config_path = os.path.abspath(args.c)
threads = args.t
output = os.path.abspath(args.o)
memory = args.m
slurm_partition = args.s

os.makedirs(output, exist_ok=True)

with open(config_path, 'r') as f:
    configs = json.load(f)

mounts = process_configs(configs)

docker_config_file = os.path.join(os.path.dirname(output), 'docker_config.json')

with open(docker_config_file, 'w') as f:
    json.dump(configs, f, indent=4)

volume_string = ' '.join([f'-B {v}:{k}:ro' for k, v in mounts.items()])

internal_command = f"java -Xmx{memory}g -jar /srv/TF-Prioritizer.jar -t {threads} -o /srv/wd -c /srv/input/configs.json"

external_command = f"apptainer run " \
                   f" -B {docker_config_file}:/srv/input/configs.json:ro {volume_string} -B {output}:/srv/wd:rw " \
                   f"docker://{image} {internal_command}"

slurm_prefix = f"srun -p {slurm_partition} --cpus-per-task={threads} --mem={memory}G --output={output}/slurm-%j.out " \
               f"--error={output}/slurm-%j.err" if slurm_partition else ""

if not args.i:
    os.system("docker image pull " + image)

os.system(slurm_prefix + external_command)
