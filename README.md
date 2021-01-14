# R2easyR

# Facilitates Mapping Experimental RNA Reactivity Data to Secondary Structures Drawn with R2R

## Why use R2easyR?

### Manually drawing experimental RNA reactivity data as colors on a secondary structure is time consuming and error prone. R2easyR, when used with the drawing program R2R, makes processing reactivity data, and drawing the reactivity on a secondary structure much easier. With a small amount of workup in a drawing program, the output of the mostly automated R2easyR/R2R pipline produces publication-ready figures rapidly.
 
## Gallery of examples

![R2easyR gallery](https://user-images.githubusercontent.com/63312483/96524425-3ae48e80-1246-11eb-9240-94b50aa1b3ec.png)

### The R2easyR package allows you to programatically convert experimental reactivity data into input files for R2R. R2R (https://sourceforge.net/projects/weinberg-r2r/) is fantastic software for drawing RNA secondary structures. R2R structures are uniform, making secondary structures drawn with R2R intuitive and publication ready. Moreover, R2R excels at displaying conservation information for RNA motifs, the primary purpose of R2R. Unfortunately, R2R contains fewer options for drawing experimental reactivity data on RNA secondary structures. While all the code R2R needs to map continuous reactivity data to colors on a secondary structure is present, R2R does not have good programatic options for supplying drawing instructions for reactivity data. Thus, R2easyR is designed to link reactivity data and the R2R syntax. R2easyR is a toolbox of R commands that integrate common RNA-lab data formats with appropriate color palettes, and write files that R2R uses to color nucleotides on secondary structuree. This approach is trivial for short RNA (<100N), and reasonably achievable for long RNA (100–3000 N).

## Overview

### The core of the R2easyR is four functions: “r2easyR.pallets”, “r2easyR.color”, “r2easyR.write”, and "r2easyR.stem_editor". “r2easy.palettes” is a palette generation function that makes 59 palettes that can be used to color reactivity data on RNA secondary structures. “r2easyR.color” is an assignment function that assigns colors to nucleotides based on the reactivity at that nucleotide. If the nucleotide reactivity is higher, the color “r2easyR.color” assigns will be more intense. “r2easyR.write” is a file writing function that writes two files (a Stockholm file and a R2R meta file) that tell R2R how to color each nucleotide on a secondary structure. "r2easyR.stem_editor" is a R2R Stockholm file editor, that supplies instructions to a stockholm file that optimizes the layout of stems in a multistem structure so that stems are less likely to clash.


![Figure_1](https://user-images.githubusercontent.com/63312483/96525110-1d182900-1248-11eb-8b60-af187eb1fdb6.png)


### The four core R2easyR functions can be used in a R script to generate generate a single RNA secondary structure in a couple minutes starting with a .csv file (Video turotial 1) or a .ct file (Video tutorial 2). Alternatively, if you are good at R scripting, one can run the four core R2easyR on a loop to print RNA secondary structures depictions in bulk. Two such bulk structure printers are built into R2easyR. The first is called "ry.tRNA" and maps reactivity data to near-"cloverleaf" tRNA structures in bulk starting with .ct and .shape formated secondary structure/reactivity data. The second, called "r2easyR.go_fast" maps reactivity data to any RNA secondary structure in bulk starting with .ct and .shape formated secondary structure/reactivity data (Video tutorial 3).

##Video tutorials

### 1.) Using R2easyR to map reactivity data to a simple RNA 2-structure using a .csv

https://youtu.be/tLZDegQ6b20

### 2.) Using R2easyR to map reactivity data to a simple RNA 2-structure using a CT file from RNAstructue

https://youtu.be/2wCFcz3NEiw

### 3.) Using R2easyR to map reactivity data to RNA 2-structures in bulk using CT files from RNAstructure

https://youtu.be/SO_8z9tvXK4

## Installation

### Prior to installing R2easyR, please install R. If you are not used to R or other command- line programs, I strongly recommend downloading and working in RStudio. To install R2easyR on your R console, open your R terminal or RStudio and type:

```{r}
install.packages(“devtools”)
```

### This will take a minute. Install any dependent packages. “devtools” is a fantastic R package for developing, distributing, and downloading R packages. I cannot guarantee that your R2easyR installation will work without “devtools”. After “devtools” is installed, type/enter:

```{r}
devtools::install_github(“JPSieg/R2easyR”)
```

### R2easyR should be installed in a few minutes. The command prompt may ask you to update dependent packages. I recommend skipping the update by entering an empty line. If you are using a fairly recent version of R, not much has changed and the updates can take a while. To check the installation when it finishes, type/enter.

```{r}
library(R2easyR)
?r2easyR.color
```

### The help files for the function “r2easyR.color” will pop up if R2easyR was installed correctly. Note the “library” function does not require quotation marks around the package. “library” loads the R2easyR package from your R library into your R memory. Thus, you will need to use “library” again after you close and restart R. You can also call R2easyR functions explicitly using the following syntax, with no need to use “library”.

```{r}
<package>::<function>
```

### For example:

```{r}
?R2easyR::r2easyR.color
```

### R2easyR relies on 4 common packages that do not come with base R: “ggplot2”, “viridis”, “RColorBrewer”, and “R.utils”. If you are a R user, you have probably already installed them. If not, or you are not sure, please install or reinstall using:

```{r}
install.packages("package")
```

### For example:

```{r}
install.packages(“ggplot2”)
```

### Reinstalling R2easyR on top of an existing R2easyR installation can corrupt the help files. Thus, if you need to reinstall R2easyR, please remove the current installation and then reinstall R2easyR. First, restart R so you are not removing a package that is currently loaded in the memory:

```{r}
q()
```

### Then restart R and run:

```{r}
devtools::install_github("JPSieg/R2easyR")
```

## Loading data into R and Formatting it for R2easyR

### The goal of loading your data into R is to supply R2easyR the nucleotide number (N), the nucleotide sequence, the secondary structure (in dot-bracket notation (.<.>.), and the Reactivity data you want to map. R2easyR functions work on data frames, a R data format that resembles a table where each column is labeled. The order of these columns does not matter, but the columns must be labeled as exactly "N", "Nucleotide", "Dotbracket", "Reactivity".

## Loading data into R with a comma separated value (.csv) text file

### This is the easiest data loading strategy, and is best for reactivity data that comes from a quantified gel. Simply open up Excel, label the first four cells (A1:A4) with the header of Figure 2, and then enter the information for your RNA. Then, save as a Comma Separated Value (“.csv”) file. For example:

![Simple_RNA](https://user-images.githubusercontent.com/63312483/104630233-88677080-5668-11eb-9925-63932e3a3fb3.png)

### Also note that missing reactivity values are specified by leaving the cell blank. Please do not use a text string like "NA" as a place holder for missing values, as this will mess with the data when it is read into R.

### Now your “.csv” file can be read into a data frame (df) using “read.csv”:

```{r}
df <- read.csv("file_name.csv")
```

### You can check your data frame (df) using the "head" function or the "View" function.

```{r}
head(df)

  N Nucleotide Dotbracket Reactivity
1 1          G          .         NA
2 2          A          .        0.6
3 3          C          <        0.2
4 4          G          <         NA
5 5          T          <         NA
6 6          A          <        0.1

View(df)
```

