# ThreeHiggs

## Container
```
(change to the directory you downloaded this software)
cd Path/To/ThreeHiggs/Directory
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

