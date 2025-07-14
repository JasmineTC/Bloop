# ThreeHiggs
All commands are to be run inside the ThreeHiggs directory
## Building a container:
```
podman build . -t threehiggs 
podman run --mount type=bind,src=$PWD,target=/ThreeHiggs -it threehiggs /bin/bash -c "cd /ThreeHiggs/src && exec /bin/bash"
```

## ~~Installating with pip (conda)~~ DEPRECIATED 
Install in developer mode with pip (again inside the ThreeHiggs directory)
```
pip install -e .
```

## Executing the code: 
```
python3 runStages.py 
```
Find cmd line flags in src/UserInput.py

## Running unit and integration tests:
```
./UnitAndIntegrationTests.sh
```
