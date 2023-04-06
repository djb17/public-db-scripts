Extracting information from KEGG to create a network of metabolites.

Directory structure:
```
kegg/
|-- other_scripts/
|-- parseKGML.py
|-- reacs.zip
|-- hsa_gene_list.json
|-- ginfo.7z
|-- hsa_list.txt
|-- README.md
```

- `other_scripts` contains of scripts used to check the information extracted.
- `parseKGML.py` python script for exactrating metabolite information
- `reacs.zip` contains the list of reactions required when parsing KGML files. Created to reduce processing time due to querying KEGG
- `hsa_gene_list.json` organizes the gene symbols associated with each KEGG gene ids.
- `ginfo.7z` contains json file that keeps track of queried results for each KEGG gene. Created to reduce processing time due to querying KEGG.
- `hsa_list.txt` is the list of all human pathways found on KEGG.
- `README.md` is the current file