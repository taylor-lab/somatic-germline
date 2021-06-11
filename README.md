# somatic-germline

Companion code to the paper: \
**_The context-specific role of germline pathogenicity in tumorigenesis_**

Code in this repository can be used to recreate the essential parts of the main figures in the manuscript. 

### Instructions
Install required R packages:
```r
install.packages(c('knitr', 'plyr', 'binom', 'readxl', 'here', 'survival', 'survminer'))
```

Get supplementary data from manuscript and download germline data from [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001858.v1.p1).

Place these in the `data` folder, `germline_mutations.maf`, and `germline_cnvs.txt`, or manually enter the file paths in the `prerequisites.R` script in the `scripts` directory.

An HTML file containing the figures can be generated from the command-line for each main figure as such:

```shell
R -e "rmarkdown::render('scripts/figure-2.Rmd', output_file = 'figure-2.html')"
```

Install python packges using pip or conda : pandas, numpy, scipt, matplotlib, seaborn

For jupyter notebooks, download notebooks and contents of the `data` folder. Ensure that the folder structure follows the same relative structure as the GitHub repository and execute notebooks.
### Citation
- pending

### Contact
E-mail any questions to [srinipr4@mskcc.org, bandlamc@mskcc.org](mailto:bandlamc@mskcc.org?subject=[GitHub]%20somatic%20germline%20paper).
