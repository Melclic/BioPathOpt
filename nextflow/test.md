Here test extracting the sink from the iML1515
```
nextflow run main.nf -entry generate_sink --input_model input/iML1515.xml
```
Test for the enrichment of the model
```
nextflow run main.nf -entry enrich --input_model input/iML1515.xml
```
