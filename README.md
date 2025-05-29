# ThreeHiggs

## Container
```
podman build . -t threehiggs 
podman run --mount type=bind,src=$PWD,target=/ThreeHiggs -it threehiggs /bin/bash
cd /ThreeHiggs/src
python3 runBenchmark.py
```

## Installation
Install in developer mode with pip:
```
pip install -e .
```

