# D-Piper - A modified Piper diagram
## Overview
D-Piper is a free source code developed in Python 3 that allows to represent large sets of hydrochemical analyses avoiding the problem of overlapping.

## How does it work?
To make a D-Piper diagram just follow these steps:
1. **Prepare the file structure:** To do so, copy and unzip to the desired location the *D_Piper_v1.zip* file.<br>Keep the file structure as shown in the figure:
![Imagen1.png](attachment:Imagen1.png)
 * The input analytical data files must be saved in the *Data* folder. The sample files used in the article to make Figures 2, 3, 4 and 5 can be found in this folder.
 * The D-Piper diagrams will be saved in the *Graphics* folder.
 * The *pycache* folder is created automatically and it is required for the *D_Piper_v1.py* script to run successfully.
 * The file *D_Piper_v1_options.txt* contains the necessary instructions and explanations to build the desired D-Piper diagram.
 * The file *D_Piper_v1.py* is the script that you have to run in Python to build a D-Piper diagram. The script is widely commented.
 * The *Tools.py* file is an auxiliary Python code that act as a module and contains all the functions with which the *D_Piper_v1.py* file operates.
 
  Do not delete or change any of the above files or folders. If necessary, you can add other folders or files. This will not affect the program execution.
 
  The script always reads the drawing options from the *D_Piper_v1_options.txt* file. If you want to save a set of options, you must do so by renaming the file or saving it to another location. For further use, you will have to name it as the original one (*D_Piper_v1_options.txt*) or update the new name in line 65 of the *D_Piper_v1.py* file.


2. **Preparing de data file:**

  The data file must be in ASCII format and saved in the *Data* folder.<br>
  The data file can have any name.<br>
  The structure of the data file should be as follows:
      * TAB-separated ASCII file
      * It must contain ten columns:
      
![Imagen2.png](attachment:Imagen2.png)

>*Identifier*: sample identification number.<br>
>*Group*: Each sample can belong to a different group. Up to nine groups can be differentiated. Although no error will be shown, setting more groups will create overlap problems.<br>
>*Eight columns with analytical results*. The ion content should be expressed in mg/L (milligrams per liter) or ppm (parts per million).


3. **Preparing the *Options* file:**

  The *D_Piper_v1_options.txt* file must be in ASCII format, separated by semicolon (;)<br>
  It is a self-explanatory file that shows all the options of the program.<br>
  It is essential to comply with the file format to avoid execution errors.<br>

  The contents of this file are explained at the end of this document. Options that can be modified are shown in $\color{red}{red}$ , explanatory text and guidelines are shown in $\color{blue}{blue}$, and variable names in **_black_**.


![dibujo.svg](attachment:dibujo.svg)
    

## Authors
All authors belong to IGME Geological Survey of Spain. C/Ríos Rosas 23, 28003 Madrid, Spain

* M. González-Jiménez         miguel.gonzalez@igme.es

* L. Moreno Merino            l.moreno@igme.es

* H. Aguilera Alonso          h.aguilera@igme.es

* E. Díaz Losada              elisabeth.diaz@igme.es


## Copyright
License: GPLv3 
