# DIA_QC_report

## Parsing Skyline reports

The `parse_data` script is used to combine .tsv reports exported from one or more Skyline documents into a sqlite database which is used in all subsequent steps. `parse_data` requires a precursor level and replicate level report from each Skyline document. The templates for both reports are in the `resources` subdirectory.

<details>
  <summary>Batch and QC database schema</summary>

![alt text](https://github.com/ajmaurais/DIA_QC_report/blob/dev/resources/schema.png?raw=true)

</details>

To create a new database:

```
parse_data --projectName <document_1> <document_1_replicate_quality.tsv> <document_1_precursor_quality.tsv> 
```

To add reports from additional documents to an existing database add the `--overwriteMode=append` flag:

```
parse_data --projectName <document_2> --overwriteMode=append <document_2_replicate_quality.tsv> <document_2_precursor_quality.tsv> 
```

## Normalzation

The `normalize_db` script performs median normalization at the precursor level and [DirectLFQ](https://github.com/MannLabs/directlfq) normalization at the protein level. The normalized values are written directly to the database. The script modifies the database in place by modifying the `normalizedArea` column in the `precursors` table and the `normalizedAbundance` column in the `proteinQuants` table.

```
normalize_db data.db3
```

## QC report

```
# Generate qmd file
generate_qc_qmd --addStdProtein 'sp|P00924|ENO1_YEAST' data.db3

# Render qmd to html
quarto render qc_report.qmd --to html
```

## Batch correction report

```
# Generate rmd file
generate_batch_rmd data.db3

# Render rmd to html
Rscript -e 'rmarkdown::render("bc_report.rmd")'
```
