import argparse
import json
import os

parser = argparse.ArgumentParser()

parser.add_argument('-c', type=str)
parser.add_argument('-t', type=str)
parser.add_argument('-o', type=str)

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


config_path = os.path.abspath(args.c)
threads = args.t
output = os.path.abspath(args.o)

os.makedirs(output, exist_ok=True)

with open(config_path, 'r') as f:
    configs = json.load(f)

mounts = process_configs(configs)

docker_config_file = os.path.join(os.path.dirname(output), 'docker_config.json')

with open(docker_config_file, 'w') as f:
    json.dump(configs, f)

volume_string = ' '.join([f'-v {v}:{k}:ro' for k, v in mounts.items()])

internal_command = f"java -jar /srv/TF-Prioritizer.jar -t {threads} -o /srv/wd -c /srv/input/configs.json"

external_command = f"docker run --user='{os.getuid()}':'{os.getgid()}'" \
                   f" -v '{docker_config_file}:/srv/input/configs.json:ro' {volume_string} -v '{output}:/srv/wd:rw,Z' " \
                   f"tfprio {internal_command}"

os.system("docker build . -t tfprio")
os.system(external_command)
