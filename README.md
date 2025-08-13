# BioPathOpt
This project provides a simulation framework for synthetic biology and metabolic 
engineering, with a focus on the analysis of metabolic models for the production 
of novel compounds. It offers methods to simplify the parsing and enrichment of 
metabolic models. Merging heterologous pathways to an existing GEM, enzyme constrained
modelling and FSEOF optimization for knockout analysis. 

# Installing

Run the following code to install

```
pip install -e .
```

You can also create a new environment file as do

## Creating the cache

The package parses the following databases: [MetaNetX](https://www.metanetx.org/) 
and [Brenda](https://www.brenda-enzymes.org/) and generates a large cache to
enrich models. The later needs to be downloaded
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

# Notebooks

Here is a non-exhaustive list of the ability of all the abilities of the package.

## Enzyme constrained model

Enzyme-constrained FBA is an extension of Flux Balance Analysis that adds limits 
on enzyme capacity by incorporating kcat values and proteomics data, leading to 
more realistic predictions of metabolic fluxes.

## FSEOF Optimization 

The FSEOF (Flux Scanning based on Enforced Objective Flux) algorithm identifies metabolic reactions whose flux increases when production of a target compound is enforced, helping to pinpoint potential overexpression targets for strain engineering.

This method find the combination of knockouts and overexpression such that there
are a few that are not missing
