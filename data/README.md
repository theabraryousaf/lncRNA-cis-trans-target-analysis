# Data folder

Put your input CSV file here.

Default expected file name:

```text
input_p0.05_lncrna_pcg.csv
```

Required format:

```text
gene_id,0h_1,0h_2,0h_3,16h_1,16h_2,16h_3,20h_1,20h_2,20h_3,24h_1,24h_2,24h_3
ENSG00000000457.14,646.58,685.77,783.36,1042.08,982.60,1038.50,1092.22,1107.30,1174.00,1553.01,1132.40,1330.04
```

The file should contain sample-wise normalized expression values, not one mean value per gene.
