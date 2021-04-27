#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#                             IGME
#                  Geological Survey of Spain
#
#                            D-Piper
#                         Version 1.0
#                          14-01-2021
#      Authors:
#		 M. González-Jiménez         miguel.gonzalez@igme.es
#		 L. Moreno Merino            l.moreno@igme.es	
#		 H. Aguilera Alonso          h.aguilera@igme.es
#		 E. Díaz Losada              elisabeth.diaz@igme.es
#		 
# 
#  If you have some problem with the code, please contact with L. Moreno Merino 
#  or M. González Jiménez

__version__    = "1.0"
__author__     = "Miguel González, Luis Moreno, Héctor Aguilera and Elisabeth Díaz"
__maintainer__ = "Miguel González and Luis Moreno"
__license__    = "GNU General Public License v3.0"

"""
-------------------------------------------------------------------------------
The purpose of this script is to draw the modified Piper diagram to represent
a large number of chemical analyses.
It reads the execution parameters located in D_Piper_v1_options.txt file and execute
all the indicated operations throughout Tools module, which must be place at the
same folder.

--> Graphic output files will be saved at 'Graphics' folder under the same name
    as the graphics title.
--> Data output files will be saved at 'Data' folder.

D_Piper_v1_options and data-entry files should be compatible with an ASCII editor
(tab separated txt files). The first one should be place in the same
folder of this script, whereas the second one should be place at 'Data' folder.

 ------------------------------------------------------------------------------
 --------------------- DATA ENTRY FILE FORMAT ---------------------------------
 ------------------------------------------------------------------------------
Position:       1         2     3  4   5   6    7     8    9   10
Variable:  Nº Analysis  Group  Ca  Mg  Na  K   HCO3  CO3  SO4  Cl

ASCII file separated by tabs   (\t)
The decimal separator is point (.)
Data file should not contain header.
-------------------------------------------------------------------------------

Hope it serves you!

"""

import Tools as Tools
from csv import reader

###############################################################################
#################------------ READING PARAMETERS ----------####################
###############################################################################
def bool_conv(x):
    ''' Convert txt boolean variables into Python real ones.
    '''
    if x == 'Y':
        return True
    elif x == 'N':
        return False
    
# The first column (Fila[0]) contains the name of the variable
# The third column (Fila[2]) contains the value of the variable
with open('D_Piper_v1_options.txt', "r") as File_Opciones:
    csv = reader(File_Opciones, delimiter=";")
    for Registro in csv:
        # General Options:
        if Registro[0] == 'Title_Piper':
            p_title = Registro[2]
        if Registro[0] == 'Title_DPiper':
            d_title = Registro[2]
        if Registro[0] == 'Bins_graph':
            bins = int(Registro[2])
        if Registro[0] == 'Point_Size':
            pt_size = int(Registro[2])
        if Registro[0] == 'Piper_D_Transparency':
            d_transparency = float(Registro[2])
        if Registro[0] == 'Piper_Transparency':
            p_transparency = float(Registro[2])
        if Registro[0] == 'Color_Unique':
            Color_Unico = Registro[2]
        if Registro[0] == 'ColorMap':
            ColorMap = Registro[2]
        if Registro[0] == 'Reverse_ColorMap':
             Reverse = bool_conv(Registro[2])

        # Group Options:
        if Registro[0] == 'GroupsNames':
            GroupsNames = Registro[2].split(',')
        # Groups colors
        if Registro[0] == 'SetColors':
            SetColores = Registro[2].split(',')
        if Registro[0] == 'Separate_Groups':
            Separate_Groups = bool_conv(Registro[2])
        # Groups to draw
        if Registro[0] == 'WhatGroupToDraw':
            WhatGroupToDraw = Registro[2].split(',')
        if Registro[0] == 'Identifier':
            Identificar = bool_conv(Registro[2])
        if Registro[0] == 'LabelColor':
            LabelColor = Registro[2]
        if Registro[0]== 'Combine':
            Combine = bool_conv(Registro[2])
        if Registro[0] == 'Combination_title':
            ov_title = str(Registro[2])

        # ---------
        # LOGARITHMIC TRANSFORMATION OF INPUT DATA
        if Registro[0] == 'Transform_Log':
            Transform_Log = bool_conv(Registro[2])

        # ---------
        # Distribution method:
        if Registro[0] == 'DivisionMethod':
            DivisionMethod = Registro[2]
        if Registro[0] == 'Intervals':
            intervals = int(Registro[2])
        # Other parameteres:
        # --------- MANUAL
        # Color scale limits are read in manual adjustment
        if Registro[0] == 'Manual_Steps_Def':
            manual_intervals = [int(i) for i in Registro[2].split(',')]
        # -------- QUANTILS
        if Registro[0] == 'InterpolationMet':
            InterpolationMet = str(Registro[2])
        # -------- STANDARD DESVIATION
        if Registro[0] == 'DSfactor':
            DStandar = float(Registro[2])
        # -------- Grid Lines ------:
        if Registro[0] == 'Grid':
            Grid = bool_conv(Registro[2])

        # ---------
        # FREQUENCY HISTOGRAM OPTIONS
        if Registro[0] == 'ShowHisto':
            ShowHisto = bool_conv(Registro[2])
        if Registro[0] == 'ColorHisto':
            ColorHisto = Registro[2]
        if Registro[0] == 'BinsHisto':
            BinsHisto = int(Registro[2])
        if Registro[0] == 'RangeHisto':
            RangeHisto = tuple(int(i) for i in (Registro[2], Registro[3]))
        if Registro[0] == 'LogHisto':
            LogHisto = bool_conv(Registro[2])

        # ---------
        # FILES OPTIONS
        if Registro[0] == 'File_Input':
            File_Input = str(Registro[2])
        if Registro[0] == 'Save_Files':
            save_file = bool_conv(Registro[2])
            
        # ---------
        # GRAPHIC FILES OPTIONS
        if Registro[0] == 'Store_Graph':
            Store_Graph = bool_conv(Registro[2])
        if Registro[0] == 'Extension_Graph':
            Extension_Graph = str(Registro[2])
        if Registro[0] == 'Graphic_Res':
            dpi = int(Registro[2])

    File_Opciones.close()  # Close file options
    
###############################################################################
##############------------ PROGRAM INITIALIZATION ------------ ################
###############################################################################
# Class Object named as 'D'
    
D = Tools.Piper(file = File_Input, sep = '\t', bins = bins)

#L ------------ Apply logarithmic transform to data if it's requiered ---------
D.logarithmic() if Transform_Log == True else None
    
###############################################################################
################------------ ASSIGNING VARIABLES ----------####################
###############################################################################
    
# Normal Piper -------------------------------------###########################
D.s_title = p_title                   # Normal Piper (s-piper) figure title.
D.s_color = Color_Unico               # Dots color for s-piper
D.pt_size = pt_size                   # Charts point size
D.s_transparency = p_transparency     # S-piper dots transparency

# D - Piper ----------------------------------------###########################
D.d_title = d_title                   # D-Piper figure title
D.color_map = ColorMap                # Color map of D-piper
D.d_transparency = d_transparency     # D-piper cells transparency
D.cmap_reversed  = Reverse            # Whether reverse the colormap or not

# Group Options ------------------------------------###########################
D.separate_groups = Separate_Groups   # Separate groups?
D.group_colors = SetColores           # Groups colors
D.group_names = GroupsNames           # Groups names       
D.groups_toDraw = WhatGroupToDraw     # Selected groups to be drawn
D.identify = Identificar           # Identify each analysis sample?
D.label_color = LabelColor            # Label color of analysis samples

# Point distribution and legend options ------------###########################
D.transform_log = Transform_Log       # Make the logarithmic transformation?
D.kind = DivisionMethod               # Distribution method   DivisionMethod
D.k = intervals                       # Distribution number of classes
D.manual_steps = manual_intervals     # User defined classes intervals
D.factor = DStandar                   # Number of standard deviations
D.q_interpolator = InterpolationMet   # Quantiles interpolation method
D.grid = Grid                         # Show grid?

# Frequency histogram ------------------------------###########################
D.show_histo = ShowHisto              # Create histogram?
D.h_color = ColorHisto                # Color of the histogram
D.h_bins = BinsHisto                  # Number of bars
D.h_range =  RangeHisto               # Tuple structure
D.h_log = LogHisto                    # Turn y-axis into log units?

# Output Graphic File Options ----------------------###########################
D.store_graph = Store_Graph           # Save output graphs?
D.dpi = dpi                           # Output graphs resolution
D.fig_format = Extension_Graph        # Output graphs format


###############################################################################
################------------ PROGRAM EXECUTION ------------ ###################
###############################################################################

#################-------------- Histogram ----------------#####################
D.histogram() if ShowHisto == True else None
#################------------ Standard Piper -------------#####################
D.s_piper_group(WhatGroupToDraw) if Separate_Groups == True else D.s_piper()
#################------------ Density Piper -------------#####################
D.method()     # Initialization of Division Method operations.

if Combine == True:
    D.d_piper()
    D.combination(ov_title, with_groups = Separate_Groups)
else:
    D.d_piper()
    
#################------------ Saving Arrays -------------######################
D.save_data() if save_file == True else None
############## --------- Relieving some RAM memory ----------- ################
del [p_title, d_title, bins,  pt_size, d_transparency, p_transparency, Color_Unico,
     ColorMap, SetColores, WhatGroupToDraw, Identificar, LabelColor, Transform_Log,
     DivisionMethod, intervals, manual_intervals, InterpolationMet, DStandar,
     Registro, ShowHisto, ColorHisto, BinsHisto, RangeHisto, LogHisto, File_Input,
     save_file, Store_Graph, Extension_Graph, dpi, Separate_Groups, GroupsNames,
     csv, File_Opciones, Combine, ov_title, Grid]

#******************************************************************************
#********************** -------- END OF PROGRAM -------- **********************
#******************************************************************************