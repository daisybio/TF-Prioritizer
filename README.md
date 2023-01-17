# Execution

The software can be executed using docker. For the following command, only [python3](https://www.python.org/downloads/),
[curl](https://curl.se/download.html) and [docker](https://docs.docker.com/get-docker/) are required.
Explanations about the configs can be found in the [config readme](configTemplates/README.md).

```bash
curl -s https://raw.githubusercontent.com/biomedbigdata/TF-Prioritizer/pipeJar/docker.py | python3 - -c [config_file] -o [output_dir] -t [threads]
```

Note, that for this approach an internet connection is required. The docker image will be downloaded from the [GitHub
Container Repository](https://ghcr.io/biomedbigdata/tfprio) on the first execution as well as with every update we
release. Furthermore, the wrapper script
will be fetched from GitHub with every execution.

If curl is not available (for example if you are using windows), or you want to be able to execute the software without
an internet connection, you can download the wrapper script
from [here](https://raw.githubusercontent.com/biomedbigdata/TF-Prioritizer/pipeJar/docker.py).

You can then execute the script using

```bash
python3 [script_path] -c [config_file] -o [output_dir] -t [threads]
```