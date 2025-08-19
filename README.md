# ThreeHiggs
**THIS CODE IS IN ACTIVE DEVELOPMENT, EXPECT BUGS AND FEATURES TO BE CHANGED AND/OR REMOVED. DOUBLE CHECK ANY RESULT FROM THE CODE.**

Download the code base with a git clone. From this point forward all commands are to be run from inside the threeHiggs directory

For easy OS compatibility we run the code inside a container from the command line. This will require podman or docker to be installed on your system. We use podman.
podman can be installed via:

- Linux/ubuntu
```sudo apt install podman```

- Windows
``` ??? ```

- Mac
```???```
  
## Using a container:
On first installisation you will need to do:
```podman build . -t threehiggs ```

With the container built we can enter the container with:
```podman run --mount type=bind,src=$PWD,target=/ThreeHiggs -it threehiggs /bin/bash -c "cd /ThreeHiggs/src && exec /bin/bash"```

This will put you in the src directory inside threeHiggs for convience. In the event of a future release using a new package you will need to rebuild the container.

From this point forward all commands are to be run from inside the threeHiggs/src directory
## Executing the code:
The code is excuted via
```python3 runStages.py ```

## Running unit and integration tests:
Ensure the code has been installed successfully via
```./UnitAndIntegrationTests.sh```

## Implementing new models:
Before we can get started with this code you need to first implement your own model in DRalgo following our examples laid out in the Mathematica directory. Store the generated txt files somewhere sensible.

Now we need to point the code to these generated text files this can be done via the command line (python3 runStages --help to see all the command line options). However, **we strongly recommend the use of a bash script or config file** for easily excuting the code in a repeatable manner. As an example we have a config file for the default Z2 case which can be found in src. 
The flags that control model dependent behaviour are: 
- Data files: --loFile etc 
- Minimisation control: --initialGuesses --var<Upper/Lower>Bounds 
- Generate benchmarks: --benchmarkFile (this should be a .py)

Generally useful flags:
- --loopOrder, - --verbose

## Generating benchmarks
We have an example benchmark generating code in src. The only thing we require from the user is the benchmark generator produces a json which we then load in benchmarkLooping.py.

