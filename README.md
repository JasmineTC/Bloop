# ThreeHiggs

## Installation
Install in developer mode with pip:
```
pip install -e .
```

## Container
A container file is also provided but it is WIP.

```
podman build . -t threehiggs 
podman run -it threehiggs /bin/bash
mv Running/runBenchmark.py src
cd Running
python3 runBenchmark.py
```

