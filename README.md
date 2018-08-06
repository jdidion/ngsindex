# ngsindex

Support for reading various NGS-related index formats:

* BAI
* TBI
* CSI


# Installation

```bash
pip install ngsindex
```

# Usage

```python
from ngsindex import IndexType, resolve_index_file, parse_index
from pathlib import Path

bam_file = Path('/path/to/reads.bam')

# Resolve the location of the BAM index file.
index_file = resolve_index_file(bam_file, IndexType.BAI)

# Load the index
index = parse_index(index_file)

# Loop through chromosome indexes
for refidx in index.ref_indexes:
    ...
```
