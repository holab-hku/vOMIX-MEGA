FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="4b7aeda06b2120349cbb37a85df1c1ccde0537f437cce501c1492d0cc5b949c5"
# Conda environment:
#   source: workflow/envs/biobakery3.yml
#   prefix: /conda-envs/6c71dbf8184f5937b984da1448b7a1cc
#   name: biobakery3
#   channels:
#     - biobakery
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - humann=3.9
#     - python=3.7.12
#     - metaphlan=4.0.6
#     - bowtie2=2.5.4
#     - diamond=2.1.10
RUN mkdir -p /conda-envs/6c71dbf8184f5937b984da1448b7a1cc
COPY workflow/envs/biobakery3.yml /conda-envs/6c71dbf8184f5937b984da1448b7a1cc/environment.yaml
# Conda environment:
#   source: workflow/envs/checkm2.yml
#   prefix: /conda-envs/8cd1f42070afc8f3535013e390ce4e5b
#   name: checkm2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - checkm2=1.0.2
#     - prodigal=2.6.3
#     - scikit-learn=0.23.2
#     - scipy=1.8.0
#     - tensorflow=2.4.0
#     - tensorflow-base=2.4.0
RUN mkdir -p /conda-envs/8cd1f42070afc8f3535013e390ce4e5b
COPY workflow/envs/checkm2.yml /conda-envs/8cd1f42070afc8f3535013e390ce4e5b/environment.yaml
# Conda environment:
#   source: workflow/envs/checkv.yml
#   prefix: /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea
#   name: checkv
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - checkv=1.0.1
#     - blast=2.15.0
#     - diamond=2.1.8
#     - hmmer=3.4
#     - python=3.10.13
RUN mkdir -p /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea
COPY workflow/envs/checkv.yml /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea/environment.yaml
# Conda environment:
#   source: workflow/envs/concoct.yml
#   prefix: /conda-envs/eb3a0523843b782c2629ca4ae6edf097
#   name: concoct
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - concoct=1.1.0
#     - python=3.12.6
RUN mkdir -p /conda-envs/eb3a0523843b782c2629ca4ae6edf097
COPY workflow/envs/concoct.yml /conda-envs/eb3a0523843b782c2629ca4ae6edf097/environment.yaml
# Conda environment:
#   source: workflow/envs/coverm.yml
#   prefix: /conda-envs/46c586ef10a6815b4b9f2830f50abbf1
#   name: coverm
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bwa=0.7.17
#     - coverm=0.6.1
#     - fastani=1.34
#     - minimap2=2.26
#     - samtools=1.9
RUN mkdir -p /conda-envs/46c586ef10a6815b4b9f2830f50abbf1
COPY workflow/envs/coverm.yml /conda-envs/46c586ef10a6815b4b9f2830f50abbf1/environment.yaml
# Conda environment:
#   source: workflow/envs/dastool.yml
#   prefix: /conda-envs/f61507bf46fad4195593a3722dbbb22b
#   name: dastool
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - das_tool=1.1.7
RUN mkdir -p /conda-envs/f61507bf46fad4195593a3722dbbb22b
COPY workflow/envs/dastool.yml /conda-envs/f61507bf46fad4195593a3722dbbb22b/environment.yaml
# Conda environment:
#   source: workflow/envs/dram.yml
#   prefix: /conda-envs/bfbcf982e543752be0db295d2d8fc7ca
#   name: dram
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.10
#     - pandas=1.5.2
#     - pytest=7.2.0
#     - scikit-bio=0.5.7
#     - prodigal=2.6.3
#     - mmseqs2==13.45111
#     - hmmer=3.3.2
#     - trnascan-se=2.0.11
#     - scipy=1.8.1
#     - sqlalchemy=1.4.46
#     - barrnap=0.9
#     - altair=4.2.0
#     - openpyxl=3.0.10
#     - networkx=2.8.8
#     - ruby=3.1.2
#     - parallel=20221122
#     - pip
#     - pip:
#       - DRAM-bio==1.5.0
RUN mkdir -p /conda-envs/bfbcf982e543752be0db295d2d8fc7ca
COPY workflow/envs/dram.yml /conda-envs/bfbcf982e543752be0db295d2d8fc7ca/environment.yaml
# Conda environment:
#   source: workflow/envs/dvf.yml
#   prefix: /conda-envs/be7220eac6b11c0146ebb95d77b0dd8b
#   name: dvf
#   channels:
#     - defaults
#   dependencies:
#     - biopython=1.78
#     - h5py=2.10.0
#     - keras=2.2.4
#     - numpy=1.19.2
#     - python=3.6.13
#     - scikit-learn=0.24.2
#     - theano=1.0.3
RUN mkdir -p /conda-envs/be7220eac6b11c0146ebb95d77b0dd8b
COPY workflow/envs/dvf.yml /conda-envs/be7220eac6b11c0146ebb95d77b0dd8b/environment.yaml
# Conda environment:
#   source: workflow/envs/eggnog-mapper.yml
#   prefix: /conda-envs/b8be70dcd988c91ea813792b098e3fd9
#   name: eggnog-mapper
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - eggnog-mapper=2.0.1
RUN mkdir -p /conda-envs/b8be70dcd988c91ea813792b098e3fd9
COPY workflow/envs/eggnog-mapper.yml /conda-envs/b8be70dcd988c91ea813792b098e3fd9/environment.yaml
# Conda environment:
#   source: workflow/envs/ete3.yml
#   prefix: /conda-envs/9410d779999acf21034e6da4819aea56
#   name: taxonomy
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - biopython=1.81
#     - ete3=3.1.3
#     - numpy=1.26.2
#     - pandas=2.1.3
#     - python=3.12.0
#     - scipy=1.11.4
RUN mkdir -p /conda-envs/9410d779999acf21034e6da4819aea56
COPY workflow/envs/ete3.yml /conda-envs/9410d779999acf21034e6da4819aea56/environment.yaml
# Conda environment:
#   source: workflow/envs/fastp.yml
#   prefix: /conda-envs/9561d60c14362ba157b9feb1fbf5df84
#   name: fastp
#   channels:
#     - anaconda
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - fastp=0.22.0
RUN mkdir -p /conda-envs/9561d60c14362ba157b9feb1fbf5df84
COPY workflow/envs/fastp.yml /conda-envs/9561d60c14362ba157b9feb1fbf5df84/environment.yaml
# Conda environment:
#   source: workflow/envs/galah.yml
#   prefix: /conda-envs/39c145009f5c11d6cc15b66771730d3c
#   name: galah
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - galah=0.4.2
RUN mkdir -p /conda-envs/39c145009f5c11d6cc15b66771730d3c
COPY workflow/envs/galah.yml /conda-envs/39c145009f5c11d6cc15b66771730d3c/environment.yaml
# Conda environment:
#   source: workflow/envs/genomad.yml
#   prefix: /conda-envs/303cc820d42025793664c0fcc95e1147
#   name: genomad
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - genomad=1.8.1
#     - keras=2.13.1
#     - tensorflow=2.13.1
RUN mkdir -p /conda-envs/303cc820d42025793664c0fcc95e1147
COPY workflow/envs/genomad.yml /conda-envs/303cc820d42025793664c0fcc95e1147/environment.yaml
# Conda environment:
#   source: workflow/envs/gtdbtk.yml
#   prefix: /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95
#   name: gtdbtk
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - gtdbtk=2.7.2
RUN mkdir -p /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95
COPY workflow/envs/gtdbtk.yml /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95/environment.yaml
# Conda environment:
#   source: workflow/envs/hostile.yml
#   prefix: /conda-envs/4ecab50c611bd20ca73e77590bda0e41
#   name: hostile
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - hostile=2.0.1
RUN mkdir -p /conda-envs/4ecab50c611bd20ca73e77590bda0e41
COPY workflow/envs/hostile.yml /conda-envs/4ecab50c611bd20ca73e77590bda0e41/environment.yaml
# Conda environment:
#   source: workflow/envs/iphop.yml
#   prefix: /conda-envs/a4751a6f51373a21124424968cc15c53
#   name: iphop
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - iphop=1.3.3
#     - python=3.8.15
RUN mkdir -p /conda-envs/a4751a6f51373a21124424968cc15c53
COPY workflow/envs/iphop.yml /conda-envs/a4751a6f51373a21124424968cc15c53/environment.yaml
# Conda environment:
#   source: workflow/envs/maxbin2.yml
#   prefix: /conda-envs/65c1faeac82fffb181455d610bd36f8b
#   name: maxbin2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - maxbin2=2.2.7
RUN mkdir -p /conda-envs/65c1faeac82fffb181455d610bd36f8b
COPY workflow/envs/maxbin2.yml /conda-envs/65c1faeac82fffb181455d610bd36f8b/environment.yaml
# Conda environment:
#   source: workflow/envs/megahit.yml
#   prefix: /conda-envs/bd6ed3d32fea9b68dbe0a22b2113f967
#   name: megahit
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - megahit=1.2.9
RUN mkdir -p /conda-envs/bd6ed3d32fea9b68dbe0a22b2113f967
COPY workflow/envs/megahit.yml /conda-envs/bd6ed3d32fea9b68dbe0a22b2113f967/environment.yaml
# Conda environment:
#   source: workflow/envs/metabat2.yml
#   prefix: /conda-envs/52f3616b8b54c3e15d0a4487cbe5315e
#   name: metabat2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - metabat2=2.17
RUN mkdir -p /conda-envs/52f3616b8b54c3e15d0a4487cbe5315e
COPY workflow/envs/metabat2.yml /conda-envs/52f3616b8b54c3e15d0a4487cbe5315e/environment.yaml
# Conda environment:
#   source: workflow/envs/metacerberus.yml
#   prefix: /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1
#   name: metacerberus
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - fastp=1.3.3
#     - fastqc=0.12.1
#     - metacerberus=1.4.0
#     - pyhmmer=0.10.15
#     - pyrodigal=3.7.1
#     - pyrodigal-gv=0.3.2
#     - python=3.12.13
#     - samtools=1.23.1
#     - scikit-learn=1.9.0
#     - scipy=1.17.1
#     - python-kaleido=0.2.1
RUN mkdir -p /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1
COPY workflow/envs/metacerberus.yml /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1/environment.yaml
# Conda environment:
#   source: workflow/envs/metaphlan.yml
#   prefix: /conda-envs/91c6a4f81e6ce3295bab6774e8fc63a4
#   name: metaphlan
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - metaphlan=4.0.6
#     - python=3.10.0
RUN mkdir -p /conda-envs/91c6a4f81e6ce3295bab6774e8fc63a4
COPY workflow/envs/metaphlan.yml /conda-envs/91c6a4f81e6ce3295bab6774e8fc63a4/environment.yaml
# Conda environment:
#   source: workflow/envs/multiqc.yml
#   prefix: /conda-envs/2627710133394bdecf097784123a7886
#   name: multiqc
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - multiqc=1.17
RUN mkdir -p /conda-envs/2627710133394bdecf097784123a7886
COPY workflow/envs/multiqc.yml /conda-envs/2627710133394bdecf097784123a7886/environment.yaml
# Conda environment:
#   source: workflow/envs/parallel-virfinder.yml
#   prefix: /conda-envs/68598a542adbb81e2995dba2c162c690
#   name: parallel-virfinder
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - parallel-virfinder=0.3.1
RUN mkdir -p /conda-envs/68598a542adbb81e2995dba2c162c690
COPY workflow/envs/parallel-virfinder.yml /conda-envs/68598a542adbb81e2995dba2c162c690/environment.yaml
# Conda environment:
#   source: workflow/envs/phabox2.yml
#   prefix: /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1
#   name: phabox2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - phabox=2.1.5
RUN mkdir -p /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1
COPY workflow/envs/phabox2.yml /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1/environment.yaml
# Conda environment:
#   source: workflow/envs/pharokka.yml
#   prefix: /conda-envs/70f56d9258420de64388e2b8d1e12709
#   name: pharokka
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - pharokka=1.9.1
#     - python=3.11.15
RUN mkdir -p /conda-envs/70f56d9258420de64388e2b8d1e12709
COPY workflow/envs/pharokka.yml /conda-envs/70f56d9258420de64388e2b8d1e12709/environment.yaml
# Conda environment:
#   source: workflow/envs/prodigal-gv.yml
#   prefix: /conda-envs/3718f9432b895f1c61eb898c6c547322
#   name: prodigal-gv
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - prodigal-gv=2.11.0
RUN mkdir -p /conda-envs/3718f9432b895f1c61eb898c6c547322
COPY workflow/envs/prodigal-gv.yml /conda-envs/3718f9432b895f1c61eb898c6c547322/environment.yaml
# Conda environment:
#   source: workflow/envs/pyhmmer.yml
#   prefix: /conda-envs/89567c28c2aec7a24f5e28aef0d15ae0
#   name: pyhmmer
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - pyhmmer=0.10.4
RUN mkdir -p /conda-envs/89567c28c2aec7a24f5e28aef0d15ae0
COPY workflow/envs/pyhmmer.yml /conda-envs/89567c28c2aec7a24f5e28aef0d15ae0/environment.yaml
# Conda environment:
#   source: workflow/envs/seqkit-biopython.yml
#   prefix: /conda-envs/b6504b6e5479ed534074060739fdf57c
#   name: seqkit
#   channels:
#     - anaconda
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - biopython=1.78
#     - numpy=1.26.0
#     - pandas=2.1.1
#     - python=3.12.0
#     - seqkit=2.6.1
RUN mkdir -p /conda-envs/b6504b6e5479ed534074060739fdf57c
COPY workflow/envs/seqkit-biopython.yml /conda-envs/b6504b6e5479ed534074060739fdf57c/environment.yaml
# Conda environment:
#   source: workflow/envs/spades.yml
#   prefix: /conda-envs/c23d8196afd2dfcd2e622801e57b851b
#   name: spades
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - spades=4.0.0
RUN mkdir -p /conda-envs/c23d8196afd2dfcd2e622801e57b851b
COPY workflow/envs/spades.yml /conda-envs/c23d8196afd2dfcd2e622801e57b851b/environment.yaml
# Conda environment:
#   source: workflow/envs/sratools-pigz.yml
#   prefix: /conda-envs/abe337f4ef1b54b5d99013c34579832e
#   name: sra-tools
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - pigz=2.8
#     - sra-tools=3.2.1
RUN mkdir -p /conda-envs/abe337f4ef1b54b5d99013c34579832e
COPY workflow/envs/sratools-pigz.yml /conda-envs/abe337f4ef1b54b5d99013c34579832e/environment.yaml
# Conda environment:
#   source: workflow/envs/strobealign.yml
#   prefix: /conda-envs/12b77e3b8e18f601c6cd28cf83cda940
#   name: strobealign
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - strobealign=0.13.0
#     - samtools=1.23.1
RUN mkdir -p /conda-envs/12b77e3b8e18f601c6cd28cf83cda940
COPY workflow/envs/strobealign.yml /conda-envs/12b77e3b8e18f601c6cd28cf83cda940/environment.yaml
# Conda environment:
#   source: workflow/envs/vamb.yml
#   prefix: /conda-envs/df6510d0002a1629524e0f02a0c8f04b
#   name: vamb
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - vamb=3.0.2
RUN mkdir -p /conda-envs/df6510d0002a1629524e0f02a0c8f04b
COPY workflow/envs/vamb.yml /conda-envs/df6510d0002a1629524e0f02a0c8f04b/environment.yaml
# Conda environment:
#   source: workflow/envs/vibrant.yml
#   prefix: /conda-envs/562ec28b54b50c320491196954064d79
#   name: vibrant
#   channels:
#       - conda-forge
#       - bioconda
#       - defaults
#   dependencies:
#       - biopython=1.79
#       - hmmer=3.4
#       - matplotlib=3.3.2
#       - numpy=1.19.5
#       - pandas=0.25.3
#       - prodigal=2.6.3
#       - python=3.7.12
#       - scikit-learn=0.21.3
#       - scipy=1.7.3
#       - seaborn=0.12.2
#       - tk=8.6.13
#       - vibrant=1.2.1
RUN mkdir -p /conda-envs/562ec28b54b50c320491196954064d79
COPY workflow/envs/vibrant.yml /conda-envs/562ec28b54b50c320491196954064d79/environment.yaml
# Conda environment:
#   source: workflow/envs/virsorter2.yml
#   prefix: /conda-envs/368c0dfb9ae7365d73597249b1529bba
#   name: virsorter2
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - virsorter=2.2.4
RUN mkdir -p /conda-envs/368c0dfb9ae7365d73597249b1529bba
COPY workflow/envs/virsorter2.yml /conda-envs/368c0dfb9ae7365d73597249b1529bba/environment.yaml

RUN conda env create --prefix /conda-envs/6c71dbf8184f5937b984da1448b7a1cc --file /conda-envs/6c71dbf8184f5937b984da1448b7a1cc/environment.yaml && \
    conda env create --prefix /conda-envs/8cd1f42070afc8f3535013e390ce4e5b --file /conda-envs/8cd1f42070afc8f3535013e390ce4e5b/environment.yaml && \
    conda env create --prefix /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea --file /conda-envs/c05eb9865829d1fa8c0fe1d9af4d4bea/environment.yaml && \
    conda env create --prefix /conda-envs/eb3a0523843b782c2629ca4ae6edf097 --file /conda-envs/eb3a0523843b782c2629ca4ae6edf097/environment.yaml && \
    conda env create --prefix /conda-envs/46c586ef10a6815b4b9f2830f50abbf1 --file /conda-envs/46c586ef10a6815b4b9f2830f50abbf1/environment.yaml && \
    conda env create --prefix /conda-envs/f61507bf46fad4195593a3722dbbb22b --file /conda-envs/f61507bf46fad4195593a3722dbbb22b/environment.yaml && \
    conda env create --prefix /conda-envs/bfbcf982e543752be0db295d2d8fc7ca --file /conda-envs/bfbcf982e543752be0db295d2d8fc7ca/environment.yaml && \
    conda env create --prefix /conda-envs/be7220eac6b11c0146ebb95d77b0dd8b --file /conda-envs/be7220eac6b11c0146ebb95d77b0dd8b/environment.yaml && \
    conda env create --prefix /conda-envs/b8be70dcd988c91ea813792b098e3fd9 --file /conda-envs/b8be70dcd988c91ea813792b098e3fd9/environment.yaml && \
    conda env create --prefix /conda-envs/9410d779999acf21034e6da4819aea56 --file /conda-envs/9410d779999acf21034e6da4819aea56/environment.yaml && \
    conda env create --prefix /conda-envs/9561d60c14362ba157b9feb1fbf5df84 --file /conda-envs/9561d60c14362ba157b9feb1fbf5df84/environment.yaml && \
    conda env create --prefix /conda-envs/39c145009f5c11d6cc15b66771730d3c --file /conda-envs/39c145009f5c11d6cc15b66771730d3c/environment.yaml && \
    conda env create --prefix /conda-envs/303cc820d42025793664c0fcc95e1147 --file /conda-envs/303cc820d42025793664c0fcc95e1147/environment.yaml && \
    conda env create --prefix /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95 --file /conda-envs/3094ed1ebef010ad0fae5d98c40a1a95/environment.yaml && \
    conda env create --prefix /conda-envs/4ecab50c611bd20ca73e77590bda0e41 --file /conda-envs/4ecab50c611bd20ca73e77590bda0e41/environment.yaml && \
    conda env create --prefix /conda-envs/a4751a6f51373a21124424968cc15c53 --file /conda-envs/a4751a6f51373a21124424968cc15c53/environment.yaml && \
    conda env create --prefix /conda-envs/65c1faeac82fffb181455d610bd36f8b --file /conda-envs/65c1faeac82fffb181455d610bd36f8b/environment.yaml && \
    conda env create --prefix /conda-envs/bd6ed3d32fea9b68dbe0a22b2113f967 --file /conda-envs/bd6ed3d32fea9b68dbe0a22b2113f967/environment.yaml && \
    conda env create --prefix /conda-envs/52f3616b8b54c3e15d0a4487cbe5315e --file /conda-envs/52f3616b8b54c3e15d0a4487cbe5315e/environment.yaml && \
    conda env create --prefix /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1 --file /conda-envs/5b969d8fbdf9ab9c8554319f85c9f1e1/environment.yaml && \
    conda env create --prefix /conda-envs/91c6a4f81e6ce3295bab6774e8fc63a4 --file /conda-envs/91c6a4f81e6ce3295bab6774e8fc63a4/environment.yaml && \
    conda env create --prefix /conda-envs/2627710133394bdecf097784123a7886 --file /conda-envs/2627710133394bdecf097784123a7886/environment.yaml && \
    conda env create --prefix /conda-envs/68598a542adbb81e2995dba2c162c690 --file /conda-envs/68598a542adbb81e2995dba2c162c690/environment.yaml && \
    conda env create --prefix /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1 --file /conda-envs/d4d369c29c6d48bdb1c3ed652ca3aed1/environment.yaml && \
    conda env create --prefix /conda-envs/70f56d9258420de64388e2b8d1e12709 --file /conda-envs/70f56d9258420de64388e2b8d1e12709/environment.yaml && \
    conda env create --prefix /conda-envs/3718f9432b895f1c61eb898c6c547322 --file /conda-envs/3718f9432b895f1c61eb898c6c547322/environment.yaml && \
    conda env create --prefix /conda-envs/89567c28c2aec7a24f5e28aef0d15ae0 --file /conda-envs/89567c28c2aec7a24f5e28aef0d15ae0/environment.yaml && \
    conda env create --prefix /conda-envs/b6504b6e5479ed534074060739fdf57c --file /conda-envs/b6504b6e5479ed534074060739fdf57c/environment.yaml && \
    conda env create --prefix /conda-envs/c23d8196afd2dfcd2e622801e57b851b --file /conda-envs/c23d8196afd2dfcd2e622801e57b851b/environment.yaml && \
    conda env create --prefix /conda-envs/abe337f4ef1b54b5d99013c34579832e --file /conda-envs/abe337f4ef1b54b5d99013c34579832e/environment.yaml && \
    conda env create --prefix /conda-envs/12b77e3b8e18f601c6cd28cf83cda940 --file /conda-envs/12b77e3b8e18f601c6cd28cf83cda940/environment.yaml && \
    conda env create --prefix /conda-envs/df6510d0002a1629524e0f02a0c8f04b --file /conda-envs/df6510d0002a1629524e0f02a0c8f04b/environment.yaml && \
    conda env create --prefix /conda-envs/562ec28b54b50c320491196954064d79 --file /conda-envs/562ec28b54b50c320491196954064d79/environment.yaml && \
    conda env create --prefix /conda-envs/368c0dfb9ae7365d73597249b1529bba --file /conda-envs/368c0dfb9ae7365d73597249b1529bba/environment.yaml && \
    conda clean --all -y
