# Multi-omics course project : Pipeline for single-cell RNA-seq and ATAC-seq data integration

## Clone the repository

Open the terminal

```
cd ~/home/CourseData/Umran
git clone https://github.com/Iguaracypsousa/brain_AD.git

```

## Download data

Data is stored in zenodo: [link] . For that, you will need to set the directory to the brain_AD and create data folder. Download the datasets required for this course into the **data** folder in **brain_AD** directory. 

```
cd brain_AD
mkdir data
```

## Setting up the environment 

Open RStudio, and click on the brain_AD.RProj to initialise the project. On the console type:

```
renv::restore()
```

