# -*- coding: utf-8 -*-
#                             IGME
#                  Geological Survey of Spain
#
#                            D-Piper
#                         Version 1.0
#                          14-01-2021
"""
Module containing all classes and functions that allows create a density Piper
diagram. All its content is executed from  D_Piper_v1 script.

*******************************************************************************
*** Object-oriented programation. Classes are interconected in a multilevel ***
*****************************   inheritance way   ***************************** 
*******************************************************************************
  _'Manipulate' (class - Grandmother)  --> Read and manipulate input data file.
  _'D_methods'  (class - Mother)       --> Data distribution functions.
  _'Piper'      (class - Daughter)     --> Plot data into Piper diagram.

Data are introduced from Piper class.
*******************************************************************************
*******************************************************************************
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.collections import PatchCollection
import warnings

###############################################################################
################## --------- DATA MANIPULATION ---------- ##################### 
###############################################################################

class Manipulate (object):
    ''' Grandmother class. Manipulate input data file. '''
        
    def __init__(self, file, sep, bins):
        '''Read the file, resolve possible format errors, transform values
        into percentage equivalentes (epm) and create matrixs arrays for each
        Piper panel.
        
        File: str
            File name (txt or csv). Eg: 'Data.txt'.
            File properties: columns separated by tabs, point as decimal
            separator, without header. Columns description:
            'ID_Analysis', 'Group', 'Ca', 'Mg', 'Na', 'K', 'HCO3', 'CO3', 'SO4',
            'Cl'
        
        sep: str, default '\t'
        bins: int
             Grid size for operations.
        '''
        
        self.file = file
        self.bins = bins
        
        #######################################################################
        ##############---------Default histogram params---------###############
        #######################################################################
        self.show_histo = True
        self.h_color = '#3EA948'   # soft green
        self.h_bins = 32           # Number of bars: auto or integer
        self.h_range = (0,80)      # Tuple: x.min  x.max if x.min=0 and x.max=0 then range is auto
        self.h_log = False         # Put the Y axis in logarithmic units?
        #######################################################################
        # Read file:
        data = pd.read_csv('{}{}'.format('Data/', file), sep = sep, header = None)
        
        # Columns arrangement:
        data.columns = ['ID_Analysis', 'Group', 'Ca', 'Mg', 'Na', 'K', 'HCO3',
                        'CO3', 'SO4', 'CL']
        info     = list(data.columns[0:2])
        cations = list(data.columns[2:6]) 
        anions  = list(data.columns[6:])
        elements = list(data.columns)[2:]
        
        # Check out of possible existing problems:
        # Comma possible problems
        data = data.applymap(lambda x: x.replace(',','.') if type(x) == str else x)
        data[elements] = data[elements].astype(float) # Data columns as integers
        data[info] = data[info].astype(str)           # Info columns as str
        
        #######################################################################
        # ----------- Percentage equivalents (epm) calculation ---------------
        #######################################################################
        
        # equivalent number (eq) = substance weight () * (valence (V) / molecular mass(Mm))
        # Conversion Factor (CF) = V/Mm ----> One per molecule
        
        Mm = {'Ca': 40.078, 'Mg': 24.305, 'Na':22.989769, 'K': 39.0983,
              'HCO3': 61.01684, 'CO3': 60.0089, 'SO4': 96.0626, 'CL':35.453}
        
        V = {'Ca': 2, 'Mg': 2, 'Na':1, 'K': 1, 'HCO3': 1, 'CO3': 2, 'SO4': 2, 'CL':1}
        
        cf = {i: V[i]/Mm[i] for i in V}
        
        # --------------------------------------------------------------------
        # --------- Convert mass units into equivalent units (eq) ------------
        # --------------------------------------------------------------------
        
        # CATIONS
        data['Ca_eq'] = (data.Ca * cf['Ca']).round(3)
        data['Mg_eq'] = (data.Mg * cf['Mg']).round(3)
        data['Na_eq'] = (data.Na * cf['Na']).round(3)
        data['K_eq']  = (data.K  * cf['K'] ).round(3)

        # ANIONS:
        data['HCO3_eq'] = (data.HCO3 * cf['HCO3']).round(3)
        data['CO3_eq']  = (data.CO3  * cf['CO3'] ).round(3)
        data['SO4_eq']  = (data.SO4  * cf['SO4'] ).round(3)
        data['CL_eq']   = (data.CL   * cf['CL']  ).round(3)
        
        [data.pop(i) for i in elements]    # Deleting useless fields
        cations = list(data.columns[2:6]) # --\
        anions  = list(data.columns[6:])  # --->  Redefined field groups
        elements = list(data.columns)[2:]  # --/
        
        # --------------------------------------------------------------------
        # --------- Convert equivalents into percentage units (epm) ----------
        # --------------------------------------------------------------------
        
        # Summatory:     
        data['CAT_sum'] = data[cations].sum(axis=1).round(3)
        data['ANI_sum'] = data[anions].sum(axis=1).round(3)
        
        # CATIONS
        data['Ca_epm']  = ((data.Ca_eq/data.CAT_sum)*100).round(3)
        data['Mg_epm']  = ((data.Mg_eq/data.CAT_sum)*100).round(3)
        data['NaK_epm'] = (((data.Na_eq + data.K_eq)/data.CAT_sum)*100).round(3)
        
        # ANIONS:
        data['HCO3CO3_epm'] = (((data.HCO3_eq + data.CO3_eq)/data.ANI_sum)*100).round(3)
        data['SO4_epm']      = ((data.SO4_eq/data.ANI_sum)*100).round(3)
        data['CL_epm']       = ((data.CL_eq/data.ANI_sum)*100).round(3)
        
        [data.pop(i) for i in elements]
        
        self.data = data
        data = self.data
        
        #-----------Creating matrix panels (Cat, An and Diamond) --------------
        Manipulate.grid(self)
                     
    def matrix(self, which):
        '''
        Creates a n x n frequency matrix, being n the size (bins) of the
        quadratic matrix.
        
        which: str. Options: ['Cation', 'Anion', 'Diamond']
          When is set to 'Cation' or 'Anion', only the lower-triangle part of
          the matrix is considered.
        '''
        
        data = self.data
        bins = self.bins
        intervals = np.linspace(0, 100, bins + 1)
        
        if which == 'Cation':
            x = list(data.NaK_epm)
            y = list(data.Mg_epm)
            
            # Creating histogram
            binned = np.histogram2d(x, y, bins = [intervals, intervals])[0]
            
            # Frequency value matrix.
            matrix = np.flip(binned, axis = 0).astype(int)
            # Generate a boolean matrix (same shape than 'matrix') and select
            # lower triangle values:
            condition = np.tril(np.ones((matrix.shape))).astype(np.bool)
            matrix    = np.where(condition, matrix, np.nan)
            matrix    = np.flip(np.transpose(matrix))

        elif which == 'Anion':
            x = list(data.CL_epm)
            y = list(data.SO4_epm)

            binned = np.histogram2d(x, y, bins = [intervals, intervals])[0]

            matrix = np.flip(binned, axis = 0).astype(int)

            condition = np.tril(np.ones((matrix.shape))).astype(np.bool)
            matrix    = np.where(condition, matrix, np.nan)
            matrix    = np.flip(np.transpose(matrix))

        elif which == 'Diamond':
            x = list(data.NaK_epm)
            y = list((100 - data.HCO3CO3_epm).round(3))
            
            binned = np.histogram2d(x, y, bins = [intervals, intervals])[0]

            matrix = np.flip(binned, axis = 0).astype(int)
            matrix = np.flip(np.transpose(matrix))

        return matrix
    
    def grid (self):
        ''' Creation of the three matrixs (by calling 'matrix' function) and the 
        resulted concatenation one (using in 'D_Methods' class). '''
        
        self.cation_mx  = Manipulate.matrix(self, which = 'Cation')
        self.anion_mx   = Manipulate.matrix(self, which = 'Anion')
        self.diamond_mx = Manipulate.matrix(self, which = 'Diamond')
        
        self.matriz = np.concatenate(
            (self.cation_mx, self.anion_mx, self.diamond_mx))
              
    def histogram(self):
        '''Create a data histogram from the concatenation matrix (self.matriz)
           It allows deciding which representation method to choose.
           Resulted figure is always saved inside 'Graphics' folder. '''
         
        h_low, h_high = self.h_range[0], self.h_range[1]
        # ---------------------- Preparing data -------------------------------
        values = self.matriz.flatten()
        values = values[~np.isnan(values)] # Remove NaN values for avoiding errors.
        values = values[values>0]          # Remove Zeros
        values = values/1000 if self.transform_log == True else values
              
        # --------------------- Histogram creation ----------------------------
        fig_h, ax_h = plt.subplots(1, figsize = (13,7))
        ax_h.hist(values, bins = self.h_bins, range = self.h_range, edgecolor='black',
                  linewidth=0.5, color = self.h_color, log = self.h_log)
        
        # ------------------- Histogram configuration -------------------------
        space = np.ceil((h_high - h_low)/10)
        x_ticks = np.arange(h_low,h_high + 1, space)
        xticklabels = (x_ticks).round(2)
        ax_h.set_xlim(h_low, h_high)
        
        ax_h.set_title('Frec. Histogram ({})'.format(self.file[0:-4]), fontsize = 26,
                       color = 'darkblue', pad = 20)
        
        if self.transform_log == True:
            ax_h.set_xlabel('Log(point density)',fontsize = 22, labelpad = 10)
            ax_h.set(xticks = x_ticks, xticklabels = xticklabels.tolist())
        else:
            ax_h.set_xlabel('point density',     fontsize = 22, labelpad = 10)
            ax_h.set(xticks = x_ticks, xticklabels = xticklabels.astype(int).tolist())
            
        ax_h.set_ylabel('Frequency', fontsize = 22, labelpad = 10)
        ax_h.tick_params(labelsize=15, pad=10, width=2, length = 10)
        ax_h.grid()

        # Make axis lines thicker:
        [ax_h.spines[axis].set_linewidth(2) for axis in ['top','bottom','left','right']]
        # ----------------------- Saving histogram ----------------------------
        fig_h.savefig('Graphics/Histogram_{}.{}'.format(self.file[:-4], self.fig_format),
                      bbox_inches = 'tight', pad_inches = 2)
        
    def logarithmic(self):
        '''
        Return log10 (x + 1) of all matrixs. Since log(1) is zero, log (x + 1)
        is computed not to lose the cells with a density of one.        
        Only computed if 'LogHisto' (in D_Piper_v1_options.txt) is set to True.
        '''
        def conversion(matriz):
            '''Logarithmic conversion. All values are multiplicated x 1000 for
            a better computing, although is not shown in D-Piper legend.'''
            # Non nan-values (First Condition)
            cond_1 = ~np.isnan(matriz) 
            # Non zero-values (Second Condition)
            cond_2 = np.ma.make_mask(matriz, copy=True, shrink=False)
            condition = np.logical_and(cond_1, cond_2)
            matriz = 1000 * np.log10(matriz, out = matriz, where = condition)
            return matriz.round(0)
        
        # Conversion to log data of all the matrixes: 
        self.matriz     = conversion(self.matriz + 1)
        self.cation_mx  = conversion(self.cation_mx + 1)
        self.anion_mx   = conversion(self.anion_mx + 1)
        self.diamond_mx = conversion(self.diamond_mx.astype(float) + 1)
        
    def save_data(self):
        ''' Save 'data' DataFrame as csv file with all columns but those
        refering to panels coordinates. Store folder: /Data'''
        
        self.data.to_csv('Data/Info_Table_{}.txt'.format(self.file[:-4]),
             index=False, sep = '\t', columns = ['ID_Analysis', 'Group', 'CAT_sum',
             'ANI_sum', 'Ca_epm','Mg_epm','NaK_epm', 'HCO3CO3_epm', 'SO4_epm',
             'CL_epm'])       
        
###############################################################################
################## --------- DIVISION METHODS ---------- ###################### 
###############################################################################
        
class D_methods(Manipulate):
    '''
     Mother class. Derived from 'Manipulate' class.
     It contains six standard methods for data classification:
     -Geometric, Manual, Equidistant, Quantiles, Standard Desviation and Jenks-
     
     There is a function per method, returning all of them an array with
     intervals limits (bp) and a str variable (text) containing layout info.
     
     'method' function acts as pointer.
     
     Class nomenclature:
        K = highest value,   L = Lowest value,  k = classes number
        bp = 'break points' of the class invertals
    '''
    
    def __init__(self, file, sep, bins):
        ''' Instantializate Grandmother class and define default parameters. '''
        Manipulate.__init__(self, file, sep, bins)
        
        #######################################################################
        ######## -------- Default Division Method parameters -------- #########
        #######################################################################
        self.kind   = 'Geometric'       # Distribution method:
        self.k      = 6                 # Number of classes
        self.factor = 0.5               # Factor for DStandard method
        self.q_interpolator = 'nearest' # Quantile interpolation method
                   
    def geometric(self):
        '''
        Geometric interval: this classification method is used for
        visualizing continuous data where distributions are skewed and
        show very pronounced rates of change (Conolly and Lake, 2006).
        
        Conolly J., Lake M. 2006. Geographical Information Systems in Archaeology
            Cambridge: Cambridge University Press.
        '''
        H, L, k = self.H, self.L, self.k
        
        r = (L/H)**(1/k)     # r: geometric factor for defining intervals
        bp = [L]
        for a in range(1, k, 1):
            limit = int(((r)**a)*int(H))
            bp.append(limit)
        bp.append(H)
        bp.sort()
        
        text = 'Geometric - Factor: {} Steps: {}'.format(round(r,2), k)
        
        return bp, text
        
    def manual(self):
        '''
        Intervals are defined by the user. This function only format input list
        and set layout information.
        A warning is shown if manual steps overpass maximum data value.
        '''
        bp = np.array(self.manual_steps).astype(int)
        text = 'Manual - Limits: {}'.format(len(bp) - 1)
        self.k = len(bp) - 1
        
        if bp.max() > self.H:
            warnings.warn('Some intervals are out of data range')
            
        return bp, text
    
    def equidistant(self):
        '''
        Creates equal intervals. This method is useful when the histogram has a
        rectangular shape (i.e. uniform distributions). The range of the data
        is divided by the number of desired classes. 
        '''
        bp = np.linspace(self.L, self.H, self.k + 1).astype(int)
        
        text = 'Equidistant - Stretches: {}'.format(self.k)
        
        return bp, text
    
    def quantiles(self):
        '''
        Creates intervals based on quantiles. This method is well suited to
        linearly distributed data. It produces irregular intervals assuring an
        equal number of values in each class and minimizing the importance of
        the class boundaries (Dent, 1999).
        
        Quantiles created throughout np.quantile() function.
        For more info, visit: 
            https://numpy.org/devdocs/reference/generated/numpy.quantile.html
        
        Dent B. D. 1999. Cartography. Thematic Map Design. Fifth Edition,
            London: McGraw-Hill Science 
        '''
        
        line = self.matrix
        line.sort()     # Matrix values are sorted in ascending order.
        #quantiles: 'q' parameter of np.quantile() function
        quantiles = np.linspace(0,1,self.k + 1).round(3)
        
        # Case in which first breaking point is 1.
        bp = np.quantile(line, q = quantiles, interpolation = self.q_interpolator)
                    
        text = 'Quantiles: - {} - Interpolation method: {}'.format(self.k, self.q_interpolator)
        
        return bp, text
    
    def DStandard(self):
        ''' Creates intervals based on the number of standard deviations.This
        would be the preferred method if data displays a normal or lognormal
        distribution (Conolly and Lake, 2006; Dent, 1999).
        Half of the resulted intervals are computed downwards from the mean value,
        and the other half are computed upwards. Always maximum and minimum 
        values are included.
        
        The number of standard deviations (self.factor) is set by 'DStandard'
        variable (in D_Piper_v1_options.txt). 
        
        A warning is shown if the resulted number of intervals is lesser of
        desired.
        
        Dent B. D. 1999. Cartography. Thematic Map Design. Fifth Edition,
            London: McGraw-Hill Science
        Conolly J., Lake M. 2006. Geographical Information Systems in Archaeology
            Cambridge: Cambridge University Press.
        '''
        
        L, H, k = self.L, self.H, np.ceil(self.k/2).astype(int)
        
        # ---------------- Stadistic and factor parameters --------------------
        std      = np.std(self.matrix)
        mean     = np.average(self.matrix)
        interval = std * self.factor  # Std relation (1, 1/2, 1/3, 1/4 ...) 
        # ----------------- Breaking point array creation ---------------------
        bp = [H]
        [bp.append(mean - interval * i) for i in range(1, k)]  # Mean-downwards
        [bp.append(mean + interval * i) for i in range(0, k)]  # Mean-upwards
        # --------- Sort and select breaking points inside data range ---------
        bp = np.array(bp).round(0).astype(int)
        bp.sort()
        bp = bp[(bp >= L) & (bp <= H)]  
        # --------- Fitting for including maximum and minimum values ----------
        bp = np.concatenate((np.array([L]), bp)) if bp[0]  > L else bp
        bp = np.concatenate((bp, np.array([H]))) if bp[-1] < H else bp

        # ------------ Labeling different for log transformed data ------------
        if self.transform_log == True:
            text = 'SD: {} - Average: {} - Classes: {}'.format(
                str(round(std/1000, 1)), str('{0:.1f}'.format(mean/1000)), len(bp) - 1)
        else:
            text = 'SD: {} - Average: {} - Classes: {}'.format(
                str(round(std, 1)), str('{0:.1f}'.format(mean)), len(bp) - 1)

        if len(bp) - 1 != self.k:
            warnings.warn("Not possible to create such amount of classes.")
            
        self.k = len(bp) - 1
        
        return bp, text
    
    def jenks(self):
        '''
        Intervals creation based on Jenks algorithm, from Jenkspy library
        (https://github.com/mthh/jenkspy).
        This data clustering algorithm optimizes the arrangement of a set of
        values into "natural" classes by seeking to minimize the average
        deviation from the class meanwhile maximizing the deviation from the
        means of the other groups. 
        '''
        try:
            from jenkspy import jenks_breaks
        except:
            raise ValueError("\nJenks library is not installed."
                "You can find it in:\nhttps://github.com/mthh/jenkspy"
                "\n------------- Installation --------------------"
                "\nFor Conda:   conda install -c conda-forge jenkspy"
                "\nFor Winpy:   pip install jenkspy")

        bp = jenks_breaks(self.matrix, nb_class = self.k) 
        bp = np.array(bp).round(0)             
        text = 'Jenks - {}  Classes'.format(self.k)
        
        return bp, text
    
    def method(self):
        '''
        Pointer function. Choose which function to execute.
        Return self.bp depending on which function is used for getting that
        breaking points.
        Default 'Geometric'
        
        A warning is shown if:
            _The classificaton method is not write correctly.
            _It's trying to input a different classificaton method.
            _The resulted number of intervals is lesser of desired.
        '''
        
        matrix = self.matriz
        matrix = matrix[~np.isnan(matrix)] # Remove NaN values for avoiding errors.
        self.matrix = matrix[matrix > 0].astype(int)  # Delete zeros.

        self.H, self.L = self.matrix.max(), self.matrix.min()

        if  self.kind == 'Geometric':
            bp, text = D_methods.geometric(self)
        elif self.kind == 'Equidistant':
            bp, text = D_methods.equidistant(self)
        elif self.kind == 'Manual':
            bp, text = D_methods.manual(self)
        elif self.kind == 'Jenks':
            bp, text = D_methods.jenks(self)
        elif self.kind == 'DStandard':
            bp, text = D_methods.DStandard(self)
        elif self.kind == 'Quantiles':
            bp, text = D_methods.quantiles(self)
        else:
            raise ValueError('There is not a distribution method called that name.')
        
        self.bp = bp
        # -------------- Modifying text for log transformed data --------------
        if self.transform_log == True:
            self.method_text = '{} (Log Transformed)'.format(text)
        else:
            self.method_text = text
        
        if len(self.bp) - 1 != self.k:
            warnings.warn("Not possible to create such amount of classes.")

        return self.bp

###############################################################################
##################### ------ PIPER DIAGRAM ------ #############################
###############################################################################
    
class Piper(D_methods):
    ''' Daughter class. Derived from 'D_methods' class.
        Creates a Piper diagram to be used as scatter plot or density plot.
    '''
    def __init__(self, file, sep, bins):
        ''' Instantializate Mother class, set diagram default params and create
        diagram coordinates for each sample.
        '''
        D_methods.__init__(self, file, sep, bins)
        
        #######################################################################
        ##############---------Default diagram params---------###############
        #######################################################################
        self.margin = 20           # Framework margin
        self.offset  = 22          # Gap between triangles and central diamond
        self.l_offset = 10         # Labels offset
        self.grid = True           # Show cells grid?
        self.sin = np.sin(np.pi/3) # sin(60º)
        self.cos = np.cos(np.pi/3) # cos(60º)
        self.tan = np.tan(np.pi/3) # tan(60º)
        # ----------- D_Piper_v1_options.txt default params -------------------
        # Scatter Diagram (S-Piper)
        self.s_title = 'Scatter Piper'  # S-Piper diagram title
        self.s_transparency = 0.1       # S-Piper diagram transparency
        self.pt_size = 2                # S-Piper diagram point size
        self.s_color = 'black'          # S-Piper diagram point color
        # Density Diagram (D-Piper)
        self.d_title = 'D - Piper'      # D-Piper diagram title
        self.color_map = 'Blues'        # D-Piper diagram color map
        self.cmap_reverse = False       # Whether reverse color palette or not
        self.d_transparency = 0.1       # D-Piper diagram cells transparency
        self.ov_title = 'Combination Diagram'
        # Output graphic options
        self.dpi = 300                  # Figure resolution (dpi, 'dots per inch')
        self.fig_format = 'png'         # Figure format
        self.store_graph = False        # Save graphics?
        self.identify = False           # Identify each sample in the diagram?
        self.label_color = 'red'        # Label color if self.identify == True
        
        # ----------- Create diagram coordinates for each sample --------------
        Piper.create_coordinates(self)
        
    def skeleton(self):
        '''Create basic structure of Piper diagram, including labels. '''

        offset, l_offset, margin = self.offset, self.l_offset ,self.margin 
        sin, cos, tan = self.sin, self.cos, self.tan

        #######################################################################
        ############ ----- BASIC SHAPE OF PIPER PLOT -----#####################
        #######################################################################
        # Triangles params:
        tri_params = {'facecolor':'None', 'lw':0.5, 'ec':'black'}
        # ---------------- DRAW CATION EQUILATERAL TRIANGLE -------------------
        cat_TRI = np.array([( 0,   0),
                            (50, sin * 100),
                            (100,  0)])
        cat_TRI = plt.Polygon(cat_TRI, **tri_params)
        
        # ---------------- DRAW ANION EQUILATERAL TRIANGLE --------------------
        an_TRI = np.array([(100 + offset,  0),
                           (150 + offset, sin * 100),
                           (200 + offset,  0)])
        an_TRI = plt.Polygon(an_TRI, **tri_params)
        
        # --------------- DRAW DIAMOND EQUILATERAL TRIANGLE -------------------
        # --------------- Conformed by fours point (H, J, I, G)----------------
        H = ((100 + offset) * 0.5        , np.cos(np.pi/6) * (100 + offset))   
        J = (((100 + offset) * 0.5 + 100), np.cos(np.pi/6) * (100 + offset)) 
        I = ((200 + offset) * 0.5        , np.cos(np.pi/6) * (200 + offset))
        G = ((200 + offset) * 0.5        , (np.cos(np.pi/6) * (200 + offset) - 
                                           2 * (np.cos(np.pi/6) * 100)))
        
        diamon_TRI = np.array([H,I,J,G])
        diamon_TRI = plt.Polygon(diamon_TRI, **tri_params)
        
        # ---------- Figure creation and panels structure addition ------------
        fig, ax = plt.subplots()
        ax.add_patch(cat_TRI)
        ax.add_patch(an_TRI)
        ax.add_patch(diamon_TRI)
        
        #######################################################################
        #######################------- Visual Grid ----------##################
        #######################################################################
        gl_params = {'lw': 0.2, 'color':'black'}
        # Grid lines set first leftward dipping, then rightwards and then horizontal
        # --------------------- Cation Triangle -------------------------------
        [plt.plot([espaciado, (espaciado/2) + 50], 
                  [0, (sin * (100 - espaciado))], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [plt.plot([100 - espaciado, -(espaciado/2) + 50], 
                  [0, (sin * (100 - espaciado))], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [plt.plot([cos * espaciado, cos * espaciado + (100 - espaciado)], 
                  [sin * espaciado, sin * espaciado], **gl_params)
                  for espaciado in range(20,100,20)]
        
        # --------------------- Anion Triangle --------------------------------
        [plt.plot([espaciado + 100 + offset, (espaciado/2) + 50 + 100 + offset], 
                  [0, (sin * (100 - espaciado))], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [plt.plot([100 - espaciado + 100 + offset, -(espaciado/2) + 50 + 100 + offset], 
                  [0, (sin * (100 - espaciado))], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [plt.plot([cos * espaciado + 100 + offset, cos * espaciado + 200 - espaciado + offset], 
                  [sin * espaciado, sin * espaciado],**gl_params)
                  for espaciado in range(20,100,20)]
        
        # --------------------- Diamond Panel ---------------------------------
        [plt.plot([((100 + offset) + espaciado) * cos, (((100 + offset) + espaciado) + 100) * cos], 
                  [sin * (offset + 100 - espaciado),
                   sin * (offset + 100 - espaciado) + sin * 100], **gl_params)
                  for espaciado in range(20,100,20)]
        
        [plt.plot([((100 + offset) + espaciado) * cos,(((100 + offset) + espaciado) + 100) * cos], 
                  [(sin * (100 + offset + espaciado)), sin * (offset + espaciado)],
                  **gl_params)
                  for espaciado in range(20,100,20)]
        
        ######################################################################
        #################---------- LABELLING ----------######################
        ######################################################################
        # Labels set by line equations.
        
        text_params = {'va':'center', 'ha':'center', 'style':'italic', 'fontsize':5}
        arrow_params = {'linewidth':0.2, 'head_width':3.5, 'head_length':4,
                        'facecolor':'black', 'edgecolor':'black'}
        x_array = np.linspace(0,50,9)
        
        # r - line:         y = tan(60º)*x + l_offset 
        ######################## -------- Mg ------ ###########################
        x, y = x_array[2],  tan * x_array[2] + l_offset
        plt.arrow(x, y, cos * 25, sin * 25, **arrow_params)
        x_txt, y_txt = x_array[2:4].mean(), sin * 22 + y
        plt.text(x_txt, y_txt, '$Mg^{2+}$', rotation = 60, **text_params)
        ######################## -------- SO4 + Cl ------ #####################
        x = x + 50 + (cos * offset)
        y = y + sin * 100 + (sin * offset)
        plt.arrow(x, y, cos * 25, sin * 25, **arrow_params)
        x_txt, y_txt = x_txt + 50 + (cos * offset), y + sin * 22
        plt.text(x_txt, y_txt, '$SO_{4}^{-2}$ + $Cl^-$', rotation = 60, **text_params)
        # r' - line:       y = tan(60º) * (x - 100 - offset)  l_offset 
        ###################### -------- CO3 + HCO3 ------ #####################
        x = x_array[4:5].mean() + 100 + offset
        y = tan * (x - 100 - offset) + l_offset
        plt.arrow(x, y, -cos * 25, -sin * 25, **arrow_params)
        x_txt, y_txt = 98 + offset + x_array[3], y - 3
        plt.text(x_txt, y_txt, '$CO_{3}^{-2}$ + $HCO_{3}^{-}$', rotation = 60, **text_params)
        # p - line:        y = tan(60º) * (-x + 100 + l_offset)  
        ######################### -------- Na + K ------ ######################
        x = np.linspace(50,100,9)[4]
        y = tan * (-x + 100 + (l_offset/tan)) 
        plt.arrow(x, y, cos * 25, -sin * 25, **arrow_params)
        x_txt,y_txt = 5 + np.linspace(50,100,9)[5:6].mean(), y -  10
        plt.text(x_txt, y_txt, '$Na^{+} + K^{+}$', rotation = 300, **text_params)
        # p' - line:       y = tan(60º) * (-x + 200 + offset + l_offset/tan) 
        ####################### -------- Ca + Mg ------ #######################
        x = 100 + offset/2 + np.linspace(0,50,9)[6]
        y = tan * x_array[2] + l_offset + sin * 100 + (sin * offset)
        plt.arrow(x, y, -cos * 25, sin * 25, **arrow_params)
        x_txt, y_txt = 99 + offset/2 + np.linspace(0,50,9)[6], y + sin * 19
        plt.text(x_txt, y_txt, '$Ca^{2+} + Mg^{2+}$', rotation = 300, **text_params)
        ######################### -------- SO4 ------ #########################
        x = 100 + offset + np.linspace(50,100,9)[6]
        y = tan * x_array[2] + l_offset
        plt.arrow(x, y, -cos * 25, sin * 25, **arrow_params)
        x_txt, y_txt = 149 + offset + x_array[6], sin * 18 + y
        plt.text(x_txt, y_txt, '$SO_{4}^{2-}$', rotation = 300, **text_params)
        # h:                 y = 0 
        ######################### -------- Ca ------ ##########################
        x, y = np.linspace(0,100,9)[5], -(cos * l_offset)
        plt.arrow(x, y, -25, 0, **arrow_params)
        x_txt, y_txt =  np.linspace(0,100,9)[4], y - 8
        plt.text(x_txt, y_txt, '$Ca^{2+}$', **text_params)
        ######################### -------- Cl ------ ##########################
        x, y = 100 + offset + np.linspace(0,100,9)[3], - (cos * l_offset)
        plt.arrow(x, y, 25, 0, **arrow_params)
        x_txt, y_txt = 100 + offset + np.linspace(0,100,9)[4], y - 8
        plt.text(x_txt, y_txt, '$Cl^{-}$', **text_params)
        
        #######################################################################
        ########### ---------- LAYOUT PROPERTIES ---------- ###################
        #######################################################################
        # Boundaries Box:
        self.box = np.array([(-margin, -margin),
                        (200 + offset + margin, -margin),
                        (200 + offset + margin, margin + self.sin * (offset + 200)),
                        (-margin, margin + self.sin * (offset + 200))])
        ax.add_patch(plt.Polygon(self.box, lw = 1, ec = 'black', fc = 'None'))

        self.ax, self.fig = ax, fig
        self.an_TRI, self.cat_TRI, self.diamon_TRI = an_TRI, cat_TRI, diamon_TRI

    def create_coordinates (self):
        ''' Create new fields with ternary coordinates of each panel. '''
        data, sin, cos = self.data, self.sin, self.cos
        
        ##################################################################
        ###################---- REPRESENTING COORDINATES ------############
        ##################################################################
        # Always coordinate center: (0,0).    
        
        # --------------------- Cation Triangle -------------------------------
        data['Cat_X'] = (data.NaK_epm + (data.Mg_epm/2)).round(3)  
        data['Cat_Y'] = (data.Mg_epm * sin).round(3)
        
        # --------------------- Anion Triangle --------------------------------
        data['An_X'] = ((data.CL_epm + data.SO4_epm/2) + (100 + self.offset)).round(3)
        data['An_Y'] = (data.SO4_epm * sin).round(3)
        
        # --------------------- Diamond panel ---------------------------------
        # Diamond:  Point out by equations systems based on dot-slope line equation.
        # Line equation:    y = mx + b   being m line's slope and b y-axis cross point.
        # m is always equals to tan(60º), with positive slope for cation triangle and 
        # negative for anion's one.
        m = sin/cos     # m = tan(60º) = sqrt(3)
        
        #          ---      Equations System:         ---
        #             y_cat = ( m * x_cat) + data.b_cat 
        #             y_an  = (-m * x_an)  + data.b_an 
        
        data['b_cat'] = data.Cat_Y  - (m * data.Cat_X)
        data['b_an']  = data.An_Y   + (m * data.An_X)  
        
        data['x_Diamond'] = ((data.b_an - data.b_cat)/(2*m)).round(3)
        data['y_Diamond'] = ((m * data.x_Diamond) + data.b_cat).round(3)
        
        [data.pop(i) for i in ['b_cat','b_an']]     # Erase useless fields
    
    def s_piper(self, for_overlap=None):
        ''' Samples plot in a scatter plot.
            
            for_overlap: bool, default None.
                Automatically set to True in 'combination' method.
        '''
        
        if for_overlap == None:
            Piper.skeleton(self)
        ax, data   = self.ax, self.data
        
        #######################################################################
        ################# --------- PIPER PLOT --------- ######################
        #######################################################################
        s_params = {'s':self.pt_size, 'alpha': self.s_transparency,
                    'color':self.s_color, 'lw':0, 'zorder':3}
        # ------------ Plotting and cliping outter part of pannels ------------
        ax.scatter(data.x_Diamond, data.y_Diamond, **s_params).set_clip_path(self.diamon_TRI)
        ax.scatter(data.An_X, data.An_Y, **s_params).set_clip_path(self.an_TRI)
        ax.scatter(data.Cat_X, data.Cat_Y, **s_params).set_clip_path(self.cat_TRI)
        
        if for_overlap == None:
            Piper.LayOut(self, figure = 's_piper')
        elif for_overlap == True:
            Piper.LayOut(self, figure = 'combination')

    def s_piper_group(self, group, for_overlap=None):
        ''' Same as s-Piper but with the option of differentiating groups.
        
            group: list, defaul (all)
            for_overlap: bool, default None
                Automatically set to True in 'combination' method.
        '''
        if for_overlap == None:
            Piper.skeleton(self)
            
        ax = self.ax
        data = self.data
        if len(group) == 1 and group[0] == '':
            raise ValueError('No group has been indicated.')
        # ---------------- Select specific group and plotting -----------------
        # ------------ Plotting and cliping outter part of pannels ------------
        kw = {'s':self.pt_size, 'alpha': self.s_transparency,'lw':0, 'zorder':3}
        
        names = {i:j for i,j in zip(group, self.group_names)}
        for i, j in zip(group, self.group_colors):
            chunk = data.loc[data.Group == str(i)]  # Selected groups
            ax.scatter(chunk.x_Diamond, chunk.y_Diamond, **kw, color = j,
              label = '{} ({})'.format(names[i], len(chunk))).set_clip_path(self.diamon_TRI)
            ax.scatter(chunk.An_X, chunk.An_Y, **kw, color = j).set_clip_path(self.an_TRI)
            ax.scatter(chunk.Cat_X, chunk.Cat_Y, **kw, color = j).set_clip_path(self.cat_TRI)
            
            if len(chunk) == 0:
                warnings.warn('Group "{}" (named "{}") does not exist.'
                    .format(i, names[i]))
                
        #-------------------------- Legend ------------------------------------
        legend = ax.legend(loc = 'upper right', bbox_to_anchor = (0.95, 0.83),
                  frameon = False, ncol = 1, fontsize = 6, handletextpad = 0.6)
        
        [i.set_sizes([13]) for i in legend.legendHandles] # Setting marker size 
        [i.set_alpha(1) for i in legend.legendHandles]    # Setting transparency
        
        if for_overlap == None:
            Piper.LayOut(self, figure = 's_piper')
        elif for_overlap == True:
            Piper.LayOut(self, figure = 'combination')

    def d_piper (self, for_overlap=None):
        '''This function takes 'bp' array calculated in D_methods class and
        generate a density diagram using Matplotlib's patch collections.
        
        The density diagram it's conformed by cells, depending its amount on the
        number of bins set. Three types of cells exist: triangular, rhomboidal
        and skewed-squares-like (called 'diamondids'), each of them created 
        through different functions. Triangular ones are placed in the right 
        artist of the cation and anion panels, being occupied the rest of the 
        surface by diamondid cells. The whole diamond panel is conformed by
        rombhus cells.
               
        Density values and cells coordinates are obtained and plotted in a
        different way depending on each panel:
         __ For cation and anion panels: by a complex of lines set from left
            vertex to top vertex and moving rightwards.
         __ For diamond panel: by a complex of lines set from top vertex to right
            vertex and moving leftwards.
            
        for_overlap: bool, default None.
            Automatically set to True in 'combination' method.
        '''
        Piper.skeleton(self)
        
        ax, bins = self.ax, self.bins
        cs = self.cs = 100/(self.bins)
        k, sin, cos = self.k, self.sin, self.cos
        cat, an, diamond = self.cation_mx, self.anion_mx, self.diamond_mx
                       
        def diamondid (x, y, anion = False):
            '''Creates an horizontal skewed-square conformed by four points,
            generated from left vertex and counting clock-wise.
            x,y: float (picked from self.data coordinates columns)
            anion: bool, default False
            '''
            if anion == True:
                x = x + 100 + self.offset
                
            diamondid = plt.Polygon(([(x           ,        y    ),
                                      (x + cs/2    , y + sin * cs),
                                      (x + cs*(3/2), y + sin * cs),
                                      (x + cs      ,       y     ) ]))
            return diamondid
        
        def triangle (x, y, anion = False):
            '''Creates equilateral an triangle conformed by three points begining
            from top vertex and going clock-wise.
            x,y: float (picked from self.data coordinates columns)
            anion: bool, default False
            '''
            if anion == True:
                x = x + 100 + self.offset
            
            tgle = plt.Polygon(([(x       ,      y      ),
                                 (x + cs/2, y - sin * cs),
                                 (x - cs/2, y - sin * cs)]))
            return tgle     
        
        def rhombus (x,y):
            '''Creates vertical rhombus conformed by four points, beginig from
            top vertex and counting clock-wise.
            x,y: float (picked from self.data coordinates columns)
            '''
            x = x + 100 + cos * self.offset
            y = y + sin * self.offset
            
            rhombus = plt.Polygon(([(x       ,       y          ),
                                  (x + cs/2  , y - sin * cs     ),
                                  (x         , y - (2 * sin * cs)),
                                  (x - cs/2  , y - sin * cs     )]))
            return rhombus
        
        #######################################################################
        ############## ------- Obtaining Cell Values ------- ##################
        #######################################################################
        # ---------------------------- Cation ---------------------------------
        condition = np.tril(np.ones((cat.shape)), k = -1).astype(np.bool)
        m_cat     = np.rot90(np.where(condition, cat, np.nan), k = 1, axes = (1,0))
        m_cat     = m_cat[~np.isnan(m_cat)]  # Flatten matrix and delete nan values.
        
        # ---------------------------- Anion ----------------------------------
        condition = np.tril(np.ones((an.shape)), k = -1).astype(np.bool)
        m_an      = np.rot90(np.where(condition, an, np.nan), k = 1, axes = (1,0))
        m_an      = m_an[~np.isnan(m_an)]   	# Flatten matrix and delete nan values.
         
        R_values = diamond.flatten()                 # Rhombus cells
        d_values = np.concatenate((m_cat, m_an))     # Diamondid cells
        t_values = np.concatenate((np.diagonal(cat), np.diagonal(an))) # Triangular cells
        
        #######################################################################
        ############ ------- Obtaining Cell Coordinates ------- ###############
        #######################################################################
        # Obtaining X coordinates of diamondids:
        X1 = np.linspace(0, 50, bins + 1)[0:-2]
        X  = list(X1)
        for col in range(1, bins):
           Xi = X1[0:-1] + cs
           [X.append(i) for i in Xi]
           X1 = Xi
        # Obtaining Y coordinates of diamondids:
        Y1 = np.linspace(0, sin * 100, bins + 1)[0:-2]
        Y  = list(Y1)
        for line in range(1, bins):
           [Y.append(i) for i in Y1[0:-1]]
           Y1 = Y1[0:-1]
         
        # Obtaining coordinates for the diagonal triangles:
        X_t = np.linspace(50,100,bins + 1)[0:-1]
        Y_t = np.linspace(sin * 100, 0, bins + 1)[0:-1]
        
        # Obtaining Rhombus coordinates
        X_d1 = np.linspace(0, 50, bins + 1)[0:-1]
        X_d = list(X_d1)
        for row in range(1, bins):
            Xi = X_d1 - cs/2
            [X_d.append(i) for i in Xi]
            X_d1 = Xi
        
        Y_d1 = np.linspace(sin * 100, 0, bins + 1)[0:-1] + (sin * 100)
        Y_d  = list(Y_d1)
        for col in range(1, bins):
            Y_di = Y_d1 - (sin * cs)
            [Y_d.append(i) for i in Y_di]
            Y_d1 = Y_di
        
        #######################################################################
        ################# ------- Patches lists ------- ####################
        #######################################################################
        # CATION:
        C_d_patches = [diamondid(x,y) for x,y in zip(X,Y)]
        C_t_patches = [triangle(x,y)  for x,y in zip(X_t,Y_t)]  
        # ANION:
        A_d_patches = [diamondid(x,y, anion = True) for x,y in zip(X,Y)]
        A_t_patches = [triangle (x,y, anion = True) for x,y in zip(X_t,Y_t)]
        # ROMBOUS:
        R_patches = [rhombus(x,y) for x,y in zip(X_d, Y_d)]
        
        #######################################################################
        ########### ----- Color Map and Classification Norm ----- #############
        #######################################################################
        cmap = plt.cm.get_cmap(self.color_map)
        palettes= ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu','GnBu', 'PuBu',
        'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
        if self.color_map in palettes:
            color_list = [cmap(i) for i in np.linspace(25,cmap.N,3).astype(int)]
            cmap = LinearSegmentedColormap.from_list('Custom map', color_list, N=cmap.N)
        
        cmap = cmap.reversed() if self.cmap_reversed else cmap
        cmap.set_under(color = 'white')
        # Manipulating boundaries for create closed intervals --> []
        bp_norm = list(np.concatenate([np.array([self.bp[0]]),
                                       np.array([i + 1 for i in self.bp[1:-1]]),
                                       np.array([self.bp[-1]])]))
        norm = BoundaryNorm(boundaries = bp_norm, ncolors = cmap.N, clip = False)
        
        #######################################################################
        ################# ------- Patches Creation ------- ####################
        #######################################################################
        if self.grid == True:    # Show the grid if desired
            col_params = {'match_original':False, 'lw': 0.1,
            'edgecolor':'#E5E5E5', 'cmap':cmap, 'norm':norm, 'alpha':self.d_transparency}
        else:  
            col_params = {'match_original':False, 'lw': 0,
            'edgecolor':'None', 'cmap':cmap, 'norm':norm, 'alpha':self.d_transparency}
        
        diamondids_coll = PatchCollection(C_d_patches + A_d_patches, **col_params)
        triangles_coll  = PatchCollection(C_t_patches + A_t_patches, **col_params)
        rombuos_coll    = PatchCollection(R_patches, **col_params)
        
        ###### ------- Adding Matrix Info to patches collections -------- #####
        diamondids_coll.set_array(d_values)
        ax.add_collection(diamondids_coll)
        
        triangles_coll.set_array(t_values)
        ax.add_collection(triangles_coll)
        
        rombuos_coll.set_array(R_values)
        ax.add_collection(rombuos_coll)
        
        #######################################################################
        ###################### --------- LEGEND ---------- ####################
        #######################################################################
        #-------- Set legend dimensions according to number of classes --------
        w, h = 0.02, 0.03         # width and height of rectangular symbols 
        height = h * (k + 1)      # Legend bar height
        # Create legend bar as a new axis [x, y, height, width]
        cbaxes = self.fig.add_axes([0.28, 0.85 - height - 0.095, w, height])
        
        if k > 11:
            raise ValueError('------ Cannot set more than 11 classes ---------')
            
        # ------------------ Add zero to legend  ------------------------------
        plt.text(1.55, 1/((k + 1) * 2),'0', fontsize = 4,
                 ha='center', va='center', transform=cbaxes.transAxes)
        
        #---------------------- Customizing Colorbar  -------------------------
        # diamond_coll collection compute for all values since it has related norm and cmap
        colorbar = self.fig.colorbar(diamondids_coll, cax = cbaxes, extend = 'min',
          extendrect = True, extendfrac = 'auto', spacing = 'uniform', drawedges = True)
        
        colorbar.outline.set_linewidth(0.3) # Contour colorbar lines.
        colorbar.ax.tick_params(size=0)     # Delete spines
        
        # ---------- Locate ticks in the middle of the intervals --------------
        up, down  = bp_norm, list(np.concatenate([[0], bp_norm[0:-1]]))
        middle = np.array([up, down]).mean(0)
        colorbar.set_ticks(middle)
        # ------------------------- Ticks format ------------------------------
        up = list(np.concatenate([np.array([self.bp[0]]), np.array(up[1:])]))
        ticks = []
        
        if self.transform_log == True:    # Log transfomed case
            for i,j in zip(up[0:-1], self.bp[1:]):
                if i != j:
                    points = np.sum(self.matrix[((i <= self.matrix) &
                                                 (j >= self.matrix))])/1000
                    percentage = (points/(np.sum(self.matrix)/1000)) * 100
                    ticks.append('{} - {}  ({}%)'.format(round(i/1000,1),
                                 round(j/1000,2), round(percentage, 1)))
                else:
                    points = np.sum(self.matrix[self.matrix == i])/1000
                    percentage = (points/(np.sum(self.matrix)/1000)) * 100
                    ticks.append('{}  ({}%)'.format(round(i/1000,1),
                                 (round(percentage, 1))))
        else:  # Non-log transfomed case
            for i,j in zip(up[0:-1], self.bp[1:]):
                if i != j:
                    points = np.sum(self.matrix[((i <= self.matrix) & 
                                                 (j >= self.matrix))])
                    percentage = (points/np.sum(self.matrix)) * 100
                    ticks.append('{} - {}  ({}%)'.format(int(i), int(j),
                                 round(percentage, 1)))
                else:
                    points = np.sum(self.matrix[self.matrix == i])
                    percentage = (points/np.sum(self.matrix)) * 100
                    ticks.append('{}  ({}%)'.format(int(i), (round(percentage, 1))))
        
        colorbar.set_ticklabels(ticks)          # Define intervals
        colorbar.ax.tick_params(labelsize = 4)  # Label's fontsize
        
        # --------------------- Legend Title ----------------------------------
        plt.text(0, 1.05 + 0.2/self.k, 'Samples per cell (samples percentage)',
                 fontsize = 4, transform = cbaxes.transAxes)
        
        if len(np.unique(self.bp)) != len(self.bp):
            warnings.warn('Intervals are too short. Try to define less intervals')
        
        #---------------------------- Text Info -------------------------------
        ax.text(-self.margin/2, self.margin + self.sin * (self.offset + 185),
                self.method_text, fontsize = 6)
        
        if not for_overlap == True:
            Piper.LayOut(self, figure = 'd_piper')
    
    def combination(self, ov_title, with_groups = False):
        '''
        Create a D -Piper and overlap it with a scatter plot (S-Piper)
        ov_title: str, default 'Combination Diagram'
                  Resulted diagram title.
        with_groups: bool, default None
                  Wether to differenciate by groups or not.
        '''
        if with_groups == False:    
            self.ov_title = ov_title
            Piper.d_piper(self, for_overlap = True)
            Piper.s_piper(self, for_overlap = True)
            
        if with_groups == True:
            self.ov_title = ov_title
            Piper.d_piper(self, for_overlap = True)
            Piper.s_piper_group(self, self.groups_toDraw, for_overlap = True)
            
    def LayOut (self, figure):
        '''    
        Identify samples (optional), set layout params (title and text)
        according to the output figure and save it (optional).
        
        figure: str. 
        '''
        ax, fig = self.ax, self.fig 
        data, margin = self.data, self.margin
        
        #######################################################################
        ############# -------- Identify points if desired --------- ###########
        #######################################################################

        if self.identify == True:
            mk_params = {'fontsize':5, 'color':self.label_color}
            [ax.annotate(i, (j,k), **mk_params) for i,j,k in zip(data.ID_Analysis,
             data.Cat_X, data.Cat_Y)]   # For Cation
            [ax.annotate(i, (j,k), **mk_params) for i,j,k in zip(data.ID_Analysis,
             data.An_X, data.An_Y)]   # For Anion
            [ax.annotate(i, (j,k), **mk_params) for i,j,k in zip(data.ID_Analysis,
             data.x_Diamond, data.y_Diamond)]   # For Diamond
            
        #######################################################################
        ################--------- TITLE AND TEXT----------#####################
        #######################################################################
        x, y = self.box[2:4].mean(0)[0], self.box[2:4].mean(0)[1] + 8
        
        title_params = {'fontsize':12, 'color':'darkblue', 'ha':'center'}
        
        if figure == 'd_piper':
        #------------------------------- Title --------------------------------
            ax.text(x, y, self.d_title, **title_params)
        
        elif figure == 's_piper':
        #------------------------------- Title --------------------------------
            ax.text(x, y, self.s_title, **title_params)
        #---------------------------- Text Info -------------------------------
            text = 'Total number of analyses: {}'.format(len(data))
            ax.text(-margin/2, margin + self.sin * (self.offset + 185), text, fontsize = 6)
        
        elif figure == 'combination':  
        #------------------------------- Title --------------------------------
            ax.text(x, y, self.ov_title, **title_params)
            
        #------------------- Proportionality and axis -------------------------
        ax.set_aspect('equal')
        ax.axis('off')
        plt.show()
        #------------------------- Saving Figure ------------------------------
        if self.store_graph == True:
            if figure == 'd_piper':
                fig.savefig('Graphics/{}.{}'.format(self.d_title, self.fig_format), dpi=self.dpi)
            elif figure == 's_piper':
                fig.savefig('Graphics/{}.{}'.format(self.s_title, self.fig_format), dpi=self.dpi)
            elif figure == 'combination':
                fig.savefig('Graphics/{}.{}'.format(self.ov_title, self.fig_format), dpi=self.dpi)
                
        #------------------------- Clearing info ------------------------------
        ax.clear()
        plt.close()
    
###############################################################################
###################---------- Custom Warnings ---------########################
###############################################################################
 
def format_warning(msg, *args, **kwargs):
    return '\nWarning: {}\n'.format(msg) # ignore everything except the message

warnings.formatwarning = format_warning
