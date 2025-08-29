# BioPathOpt

BioPathOpt is a Python package designed to streamline metabolic engineering workflows by integrating advanced modeling and optimization tools into a unified framework. It enables researchers to enrich stoichiometric metabolic models with additional biochemical data, seamlessly incorporate enzyme-constrained modeling approaches, and perform targeted optimization for the production of desired molecules.

## Enzyme Constrained Modelling

Enzyme-constrained modelling extends traditional metabolic models by limiting reaction rates based on enzyme capacities. 
This adds biological realism, captures proteome allocation trade-offs, and improves prediction accuracy for metabolic engineering.
There are two different modelling frameworks for this type of modelling [GECKO](https://github.com/SysBioChalmers/GECKO) written in 
MATLAB and [ECMpy](https://github.com/tibbdc/ECMpy). The former is impractical if you do not have the right license and the later
is slow and has not been updated for a while. To this end, we provide this package to help the metabolic engineering community
run these types of models with more ease.

Below are plots that illustrate the improvements of this type of modelling framework to accurately simulate cell metabolism

|![](notebooks/img/original_measured_v_predicted.png)|![](notebooks/img/ecModel_measured_v_predicted.png)|
|:-:|:-:|
|Original iML1515 model predictions for various carbon sources|EC iML1515 model predictions for various carbon sources|

We see that there are clear improvements:

![](notebooks/img/prediction_improvements.png)

# Docker

The easiest way to use this project is to use the docker functionality that will
open a notebook that enables you to run all aspects of this project. Below are 
three options:

## Running with Docker Compose

This project ships with a flexible `docker-compose.yml` that lets you run **BioPathOpt** in three different ways:

1. **Use a prebuilt image from Docker Hub**  
2. **Build locally from `images/Dockerfile.cache`**  
3. **Build locally from `images/Dockerfile` with a custom `BRENDA_FILE` argument**  

### 1. Run from Docker Hub (recommended)

Pulls the image directly from Docker Hub:

```bash
docker compose --profile hub up --build
```

### 2. Build from `images/Dockerfile.cache`

Uses the cached build recipe for faster local builds. This is because the repeated
pubchem REST requests will make the image fail. We recommend that you use the
following file to generate the cache before :

```bash
docker compose --profile cache up --build
```

### 3. Build from `images/Dockerfile` with `BRENDA_FILE`

Allows you to specify which **Brenda flatfile** to include at build time. Note
that this will more likely fail because of the repeated REST requests to pubchem
and we recommend that you either generate the cache locally, or use the image on 
dockerhub:

```bash
# Use default brenda_2024_1.json
docker compose --profile full up --build

# Or override with a custom file
BRENDA_FILE=brenda_2025_1.json docker compose --profile full up --build
```

### Notes

- The JupyterLab server will be available at **http://localhost:8085**.  
- Notebooks in your local `./notebooks/` directory are mounted inside the container at `/home/notebooks/`.  
- No authentication token or password is set by default (`--NotebookApp.token='' --NotebookApp.password=''`).  

# Installing

Run the following code to install by travelling to the github folder and running:

```
pip install -e .
```

# Creating the cache

> [!IMPORTANT]
> You must manually download the brenda JSON on their website manually

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

Here is a non-exhaustive list of the features of this package.

## Enrich model

The [Read model notebook](notebooks/read_model.ipynb) shows you how to pass any cobrapy model
and enrich the annotations automatically

## Enzyme constrained model

[Enzyme-constrained notebook](notebooks/build_enzyme_constrained_model.ipynb) shows you 
how to build an enzyme constrained model and using an E.Coli model demonstrates the 
improvements in this modelling method vs the classical model by comparing it to
actual measured data. (taken from [ECMpy](https://github.com/tibbdc/ECMpy))

## FSEOF Optimization 

The [FSEOF notebook]() (Flux Scanning based on Enforced Objective Flux) algorithm identifies 
metabolic reactions whose flux increases when production of a target compound is 
enforced, helping to pinpoint potential overexpression targets for strain engineering.

# Acknowledgements 

This package would not be possible without the work of [ECMpy](https://github.com/tibbdc/ECMpy) 
and [FVSEOF](https://github.com/LucasCoppens/fvseof/tree/main) where we I copied
and altered a lot of the code. This is built relying heavily on [CobraPy](https://github.com/opencobra/cobrapy).

# References

- Mao, Zhitao, et al. "ECMpy 2.0: A Python package for automated construction and analysis of enzyme-constrained models." Synthetic and systems biotechnology 9.3 (2024): 494-502.
- Ebrahim, Ali, et al. "COBRApy: constraints-based reconstruction and analysis for python." BMC systems biology 7.1 (2013): 74.
- Chen, Yu, et al. "Reconstruction, simulation and analysis of enzyme-constrained metabolic models using GECKO Toolbox 3.0." Nature protocols 19.3 (2024): 629-667.

