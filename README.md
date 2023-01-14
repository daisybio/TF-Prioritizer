# Execution

The software can be executed using docker. For the following command, only python3, curl and docker are required.

```{bash, warning=FALSE, message=FALSE}
curl -s https://raw.githubusercontent.com/biomedbigdata/TF-Prioritizer/pipeJar/docker.py | python3 - -c [config_file] -o [output_dir] -t [threads]
```