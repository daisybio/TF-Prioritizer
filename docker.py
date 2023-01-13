import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-c', type=str)
parser.add_argument('-t', type=str)
parser.add_argument('-o', type=str)

args = parser.parse_args()

config = os.path.abspath(args.c)
threads = args.t
output = args.o

internal_command = f"java -jar /srv/TF-Prioritizer.jar -t {threads} -o /srv/wd -c /srv/input/configs.json"

external_command = f"docker-compose run --user='{os.getuid()}':'{os.getgid()}'" \
                   f" -v '{config}:/srv/input/configs.json' tfprio {internal_command}"

os.system(external_command)
