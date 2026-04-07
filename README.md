# BLOOP (Beyond one LOOp Phase transition)
** THIS CODE HAS BEEN MOVED TO https://github.com/MOREHIGGS/Bloop - this page will remain to preserve dev history **

## Installing the code:
Download the code base with a git clone. From this point forward all commands are to be run from inside the Bloop directory

For cross platform compatibility and clean installation environment we recommend install the code in a container using podman (or docker). Alternatively, the code can be pip installed locally.
<details>
 <summary>Installition with podman</summary>

If podman is not already installed on your computer you can follow the relavent instructions:

<details>
 
<summary>Windows (using winget in the powershell)</summary>
Virtualization must be enabled in the BIOS (default for windows 11) 

 ```bash
winget install -e --id RedHat.Podman
```
```bash 
podman machine init
```
```bash
podman machine start
 ```

</details>

<details>
<summary>macOS</summary>
Virtualization must be enabled (on by default)
 
```bash 
brew install podman
```
```bash
podman machine init
```
```bash
podman machine start
```

</details>

<details>
<summary>Linux(Ubuntu/Debian)</summary>


```bash 
sudo apt install podman
```

</details>

If the container hasn't already been built then run (from inside the Bloop directory):

```bash 
podman build . -t bloop
```
The container needs to only be built once unless a new version of BLOOP requires the container to be rebuilt (this will be communicated via the release notes). With the container built we can enter the container with:

```bash 
podman run --mount type=bind,src=$PWD,target=/Bloop -it bloop /bin/bash -c "cd /Bloop/Run && exec /bin/bash"
```

This will put you in the Run directory inside Bloop for convience. 
</details>

<details>
<summary>Installing with PIP </summary>
From inside the Bloop directory run:

 ```bash
pip install -e .
```
 
</details>

## Using BLOOP:
Bloop can be run with default settings by excuting :

 ```bash
bloop
```
 
from the cmd line. Default settings for bloop run the Z2-3HDM model. Details of which can be found here https://arxiv.org/abs/2511.04636. 

BLOOP is controlled by command line arguments. To see a full list of arguments simply run  

```bash
bloop --help
```
A brief outline of the most important flags:

- Model file paths: --modelDirectory (relative to Build) --loFilePath (relative to modelDirectory)
- Generating benchmarks: --bmGeneratorModule
- --loopOrder
- --verbose
- --configFilePath (relative to cwd)
- --bPlot


## Implementing new models:
Before we can get started with this code you need to first implement your own model in DRalgo following our examples laid out in the Mathematica directory. 
Store the generated txt/JSON files inside build.

BLOOP needs to be given the model file paths via the cmd line. To make this easier **we strongly recommend the use of a bash script or config file**. We have an example config file for the default Z2 case which can be found in Run. 
  
We have an example benchmark generating code in Source. The only thing we require from the user is the benchmark generator produces a json which we then load in benchmarkLooping.py.

