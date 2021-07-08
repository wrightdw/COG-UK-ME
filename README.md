# COG-UK-Dashboard

## Pipeline
The R Markdown file mutations.Rmd forms the pipeline that takes the large CSV files containing sequence metadata, mutations and deletions from the MSA daily dataset on CLIMB and prepares the pre-computed data files that back the COG-ME dashboard.

The YAML section at the top of mutations.Rmd defines the dependencies of the pipeline and the data of the MSA dataset to be processed. To change the date, edit the params: dataset-date field in the format YYYY-MM-DD 
Save mutations.Rmd and convert it to a plain R script by running knitr::purl("mutations.Rmd", documentation = 1)
This produces file mutations.R. Transfer mutations.R to CLIMB using SFTP. 

### CLIMB

Conda is used to create a local R environment isolated from the native R installation on CLIMB in order to manage R and its dependencies.

### Install miniconda
From the command line on CLIMB
Download miniconda from https://docs.conda.io/en/latest/miniconda.html
Choose Linux installers > Latest Python version > Miniconda 3 Linux 64-bit
Download from the link on the page e.g. wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-MacOSX-x86_64.sh
Run the installer shell script

#### Create conda environment

conda env create --file environment.yml

This installs R and most of the required R packages

Package seqinr is not installed via conda as the bioconda version is out of date. 
Install seqinr manually at the R command line using install.packages

#### Pipeline dependencies
Change to the working directory
cd /cephfs/covid/bham/climb-covid19-wrightd/COG-UK-Dashboard

Upload file dependencies to the working directory
Dependencies:

From GitHub repo SARS-CoV-2-Antigenic:
spike_escape_info.csv
spike_escape_extra.csv
tcell_long.csv
Remdesivir Mutants (directory and all files below)

From /cephfs/covid/bham/results/msa/latest/metadata/
cog_global_DATE_consortium.csv
cog.deletions.tsv

From GitHub repo COG-UK-Dashboard:
Prediction_from_Morten.txt
sequence.fasta
VUI and VOC.csv

Create a subdirectory named C0G-UK if not already created

#### Running the pipeline
The pipeline is run via the shell script run_mutations.sh which may be automated using crontab.
Log file run_mutations.log
When the pipeline has completed, the output is in subdirectory COG-UK
Automatically generated output may be downloaded from directory /cephfs/covid/bham/climb-covid19-wrightd/COG-UK-Dashboard/COG-UK
Tranfer the file vui_voc.rds and the sdirectory named according to the date of the dataset e.g. 2021-06-22
to your local computer.

### Local
Edit global.R
Update the date in the line near the top:
dataset_date <- ymd("2021-06-22")
Run the Shiny app from RStudio to verify.

Deploy to web server sars2
SSH into sars2.cvr.gla.ac.uk

sudo -i
cd  /opt/shiny-server/my-apps/COG-UK/
systemctl stop shiny-server.service
mv /home/userid/COG-UK/* .
systemctl start shiny-server.service
