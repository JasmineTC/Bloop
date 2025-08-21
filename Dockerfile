FROM python:3

ENV PIP_ROOT_USER_ACTION=ignore

RUN pip install --upgrade pip

RUN pip install --no-cache-dir pdg 
RUN pip install --no-cache-dir numba 
RUN pip install --no-cache-dir nlopt 
RUN pip install --no-cache-dir numpy 
RUN pip install --no-cache-dir scipy 
RUN pip install --no-cache-dir ijson 
RUN pip install --no-cache-dir pathos 
RUN pip install --no-cache-dir sympy 
RUN pip install --no-cache-dir matplotlib 
RUN pip install --no-cache-dir importlib
RUN pip install --no-cache-dir line_profiler
