FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="ab81d9490037ebb43089560a1d9e46f99d4accc889946ed2a15e5ec84ab17f6d"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/dvf.yml
#   prefix: /conda-envs/b4951c3ba8637e2c58bba041fe8881b0
#   name: dvf
#   channels:
#     - defaults
#   dependencies:
#     - biopython=1.78=py36h7f8727e_0
#     - h5py=2.10.0=py36hd6299e0_1
#     - keras=2.2.4=0
#     - numpy=1.19.2=py36h54aff64_0
#     - python=3.6.13=h12debd9_1
#     - scikit-learn=0.24.2=py36ha9443f7_0
#     - theano=1.0.3=py36hfd86e86_0
RUN mkdir -p /conda-envs/b4951c3ba8637e2c58bba041fe8881b0
COPY workflow/envs/dvf.yml /conda-envs/b4951c3ba8637e2c58bba041fe8881b0/environment.yaml

# Conda environment:
#   source: workflow/envs/genomad.yml
#   prefix: /conda-envs/7173b31b5691e611d5cfe5bfa202934b
#   name: genomad
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - genomad=1.8.1=pyhdfd78af_0
#     - keras=2.13.1=pyhd8ed1ab_0
#     - tensorflow=2.13.1=cpu_py310hd1aba9c_1
RUN mkdir -p /conda-envs/7173b31b5691e611d5cfe5bfa202934b
COPY workflow/envs/genomad.yml /conda-envs/7173b31b5691e611d5cfe5bfa202934b/environment.yaml

# Conda environment:
#   source: workflow/envs/parallel-virfinder.yml
#   prefix: /conda-envs/d77a208a7f81cb587f2d3b0462212a92
#   name: parallel-virfinder
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - parallel-virfinder=0.3.1=py310hdfd78af_0
RUN mkdir -p /conda-envs/d77a208a7f81cb587f2d3b0462212a92
COPY workflow/envs/parallel-virfinder.yml /conda-envs/d77a208a7f81cb587f2d3b0462212a92/environment.yaml

# Conda environment:
#   source: workflow/envs/phabox2.yml
#   prefix: /conda-envs/f24788767f0fd430ce4277dd62814291
#   name: phabox2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - phabox=2.1.5=pyhdfd78af_0
RUN mkdir -p /conda-envs/f24788767f0fd430ce4277dd62814291
COPY workflow/envs/phabox2.yml /conda-envs/f24788767f0fd430ce4277dd62814291/environment.yaml

# Conda environment:
#   source: workflow/envs/seqkit-biopython.yml
#   prefix: /conda-envs/6b8a84643dc3ec104144189c2529c0a7
#   name: seqkit
#   channels:
#     - anaconda
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - biopython=1.78=py312h5eee18b_0
#     - numpy=1.26.0=py312hc5e2394_0
#     - pandas=2.1.1=py312h526ad5a_0
#     - python=3.12.0=h996f2a0_0
#     - seqkit=2.6.1=h9ee0642_0
RUN mkdir -p /conda-envs/6b8a84643dc3ec104144189c2529c0a7
COPY workflow/envs/seqkit-biopython.yml /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml

# Conda environment:
#   source: workflow/envs/vibrant.yml
#   prefix: /conda-envs/92ad042059227bff0dc45a8e9deb74fc
#   name: vibrant
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - biopython=1.79=py37h540881e_2
#       - hmmer=3.4=hb6cb901_4
#       - matplotlib=3.3.2=0
#       - numpy=1.19.5=py37h3e96413_3
#       - pandas=0.25.3=py37hb3f55d8_0
#       - prodigal=2.6.3=h577a1d6_11
#       - python=3.7.12=hf930737_100_cpython
#       - scikit-learn=0.21.3=py37hcdab131_0
#       - scipy=1.7.3=py37hf2a6cf1_0
#       - seaborn=0.12.2=hd8ed1ab_0
#       - tk=8.6.13=noxft_hd72426e_102
#       - vibrant=1.2.1=hdfd78af_4
RUN mkdir -p /conda-envs/92ad042059227bff0dc45a8e9deb74fc
COPY workflow/envs/vibrant.yml /conda-envs/92ad042059227bff0dc45a8e9deb74fc/environment.yaml

# Conda environment:
#   source: workflow/envs/virsorter2.yml
#   prefix: /conda-envs/ecca5f750a65e1a1b5e57fa7da19e99d
#   name: virsorter2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - virsorter=2.2.4=pyhdfd78af_1
RUN mkdir -p /conda-envs/ecca5f750a65e1a1b5e57fa7da19e99d
COPY workflow/envs/virsorter2.yml /conda-envs/ecca5f750a65e1a1b5e57fa7da19e99d/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/b4951c3ba8637e2c58bba041fe8881b0 --file /conda-envs/b4951c3ba8637e2c58bba041fe8881b0/environment.yaml && \
    conda env create --prefix /conda-envs/7173b31b5691e611d5cfe5bfa202934b --file /conda-envs/7173b31b5691e611d5cfe5bfa202934b/environment.yaml && \
    conda env create --prefix /conda-envs/d77a208a7f81cb587f2d3b0462212a92 --file /conda-envs/d77a208a7f81cb587f2d3b0462212a92/environment.yaml && \
    conda env create --prefix /conda-envs/f24788767f0fd430ce4277dd62814291 --file /conda-envs/f24788767f0fd430ce4277dd62814291/environment.yaml && \
    conda env create --prefix /conda-envs/6b8a84643dc3ec104144189c2529c0a7 --file /conda-envs/6b8a84643dc3ec104144189c2529c0a7/environment.yaml && \
    conda env create --prefix /conda-envs/92ad042059227bff0dc45a8e9deb74fc --file /conda-envs/92ad042059227bff0dc45a8e9deb74fc/environment.yaml && \
    conda env create --prefix /conda-envs/ecca5f750a65e1a1b5e57fa7da19e99d --file /conda-envs/ecca5f750a65e1a1b5e57fa7da19e99d/environment.yaml && \
    conda clean --all -y
