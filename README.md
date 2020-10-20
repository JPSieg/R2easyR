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


### The four core R2easyR functions can be used in a R script to generate generate a single RNA secondary structure in a couple minutes starting with a .csv file (Video turotial 1) or a .ct file (Video tutorial 2) in a couple of minutes. Alternatively, if you are good at R scripting, one can run the four core R2easyR on a loop to print RNA secondary structures depictions in bulk. Two such bulk structure printers are built into R2easyR. The first is called "ry.tRNA" and maps reactivity data to near-"cloverleaf" tRNA structures in bulk starting with .ct and .shape formated secondary structure and reactivity data. The second, called "r2easyR.go_fast" maps reactivity data to any RNA secondary structure in bulk starting with .ct and .shape formated secondary structure and reactivity data (Video tutorial 3).

## Installation

## In your R console
###    > install.packages("devtools")

###    > devtools::install_github("JPSieg/R2easyR")

## Video tutorials

### 1.) Using R2easyR to map reactivity data to a simple RNA 2-structure using a .csv

https://youtu.be/tLZDegQ6b20

### 2.) Using R2easyR to map reactivity data to a simple RNA 2-structure using a CT file from RNAstructue

https://youtu.be/2wCFcz3NEiw

### 3.) Using R2easyR to map reactivity data to RNA 2-structures in bulk using CT files from RNAstructure

https://youtu.be/SO_8z9tvXK4
