# BioPathOpt
This project provides a simulation framework for synthetic biology and metabolic 
engineering, with a focus on the analysis of metabolic models for the production 
of novel compounds. It offers methods to simplify the parsing and enrichment of 
metabolic models and the simulation, exploration of optimization of the 

# Installing

Run the following code to install

```
pip install -e .
```

You can also create a new environment file as do

## Creating the cache

The package parses the following databases: [MetaNetX](https://www.metanetx.org/) 
and [Brenda](https://www.brenda-enzymes.org/). The later needs to be downloaded
[here](https://www.brenda-enzymes.org/download.php) and compressed: 

```
from biopathopt.utils import write_compressed_json
import json
write_compressed_json(
    json.load(open('brenda_2024_1.json')), 
    'biopathopt/flatfiles/json_brenda.json.gz',
    )
```

To generate the cache run the following:

```
import biopathopt
data = biopathopt.Data()
data.refresh_cache()
```
