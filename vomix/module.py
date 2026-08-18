class Module:
    def __init__(
        self,
        module,
        workdir,
        outdir,
        datadir,
        samplelist,
        fasta="",
        fastadir="",
        sample_name="",
        assembly_ids="",
        latest_run="",
        keep_intermediates=False,
        setup_database=True,
        ncbi_email="",
        ncbi_api_key="",
        verbose=False,
        custom_config=None,
        reset=False,
    ):
        self.module = module
        self.workdir = workdir
        self.outdir = outdir
        self.datadir = datadir
        self.samplelist = samplelist
        self.fasta = fasta
        self.fastadir = fastadir
        self.sample_name = sample_name
        self.assembly_ids = assembly_ids
        self.latest_run = latest_run
        self.keep_intermediates = keep_intermediates
        self.setup_database = setup_database
        self.ncbi_email = ncbi_email
        self.ncbi_api_key = ncbi_api_key
        self.custom_config = custom_config
        self.reset = reset
