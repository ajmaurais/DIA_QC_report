# DIA_QC_report

The `dia_qc` script contains several sub-commands to generate QC and batch reports from DIA proteomics data. Reports from 1 or more Skyline documents are combined into a single sqlite database which is modified or read by subsequent commands.

## Parsing Skyline reports

The `parse` sub-command is used to combine .tsv reports exported from one or more Skyline documents into a sqlite database which is used in all subsequent steps. `parse` requires a precursor level and replicate level report from each Skyline document. The templates for both reports are in the `resources` subdirectory.

<details>
  <summary>Batch and QC database schema</summary>

![alt text](https://github.com/ajmaurais/DIA_QC_report/blob/dev/resources/schema.png?raw=true)

</details>

To create a new database:

```
dia_qc parse --projectName <document_1> <document_1_replicate_quality.tsv> <document_1_precursor_quality.tsv> 
```

To add reports from additional documents to an existing database add the `--overwriteMode=append` flag:

```
dia_qc parse --projectName <document_2> --overwriteMode=append <document_2_replicate_quality.tsv> <document_2_precursor_quality.tsv> 
```

## Normalzation

The `normalize` sub-command performs median normalization at the precursor level and [DirectLFQ](https://github.com/MannLabs/directlfq) normalization at the protein level. The normalized values are written directly to the database. The script modifies the database in place by modifying the `normalizedArea` column in the `precursors` table and the `normalizedAbundance` column in the `proteinQuants` table.

```
dia_qc normalize data.db3
```

## QC report

```
# Generate qmd file
dia_qc qc_qmd --addStdProtein 'sp|P00924|ENO1_YEAST' data.db3

# Render qmd to html
quarto render qc_report.qmd --to html
```

## Batch correction report

```
# Generate rmd file
dia_qc batch_rmd data.db3

# Render rmd to html
Rscript -e 'rmarkdown::render("bc_report.rmd")'
```

## Other sub-commands

### `db_export`

The `db_export` sub-command is used to export precursor, protein, or metadate tables from the batch/qc database to a tsv file.

```
dia_qc db_export data.db3
```

### `metadata_convert`

The `metadata_conver` sub-command is used to interconvert between all the metadata file formats which are supported by `dia_qc parse`.

```
# Convert tsv to json
dia_qc metadata_convert -o json <metadata_tsv>

# Convert tsv to Skyline annotation csv
dia_qc metadata_convert -o skyline <metadata_tsv>

# Convert Skyline annotation csv to tsv
dia_qc metadata_convert -o tsv <annotation_csv>
```
