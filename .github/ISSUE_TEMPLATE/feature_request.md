---
name: Feature request
about: Suggest a software or feature to be added to vOMIX-MEGA !
title: Feature request
labels: enhancement
assignees: ''
type: Feature

---

**Before Submitting**
- [ ] I have searched existing issues and this is not a duplicate.
- [ ] I have read the [documentation](https://vomix-mega.readthedocs.io/) to ensure this is not already implemented.
- [ ] I have verified the tool/software I am requesting is available via Conda/Bioconda.

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen. Specify the new tool, algorithm, or workflow step you propose.

**Efficiency and Benchmarking Justification**
vOMIX-MEGA strives to implement highly efficient, well-benchmarked, and widely adopted tools. Therefore, any proposed new software must be **thoroughly tested** by the requester. Requests for untested or niche software will only be considered if they provide considerable novelty or a significant performance advantage over existing tools.

Please provide any benchmark data, publications, or performance metrics that support your request.

```
[Paste any relevant benchmarking results or citations here]
```

**Conda Installability**
Since every rule in vOMIX-MEGA runs inside its own isolated Conda environment (defined in the `envs/` directory), the proposed software **must be available via Conda** (or Bioconda). Please provide the exact conda package name and version required, or a `.yml` installation file. 

```
e.g., bioconda::kraken2=2.1.3
```

**Alternatives considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context, use cases, or screenshots about the feature request here.
