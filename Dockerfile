FROM ubuntu:24.10

RUN apt update
RUN apt install -y python3
RUN apt install -y python3-numpy
RUN apt install -y python3-scipy
RUN apt install -y python3-ijson
RUN apt install -y python3-pathos
RUN apt install -y python3-numba
RUN apt install -y python3-nlopt
RUN apt install -y python3-sympy

