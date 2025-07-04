# ThreeHiggs

## Container
These commands are to be run inside the ThreeHiggs directory
```
podman build . -t threehiggs 
podman run --mount type=bind,src=$PWD,target=/ThreeHiggs -it threehiggs /bin/bash -c "cd /ThreeHiggs/src && exec /bin/bash"
python3 runBenchmark.py
```

## Installation
Install in developer mode with pip:
Depreciated(?) in favour of container
```
pip install -e .
```

