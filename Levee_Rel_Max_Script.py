#!/usr/bin/env python
# coding: utf-8

# In[1]:



#Imports
import time
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import warnings
import plotly.express as px
#import matplotlib.mlab as mlab
from scipy import optimize
import pandas as pd
from scipy import stats
from scipy.signal import argrelextrema
from scipy.stats import linregress
from scipy.optimize import curve_fit
import multiprocessing as mp
import traceback


# In[2]:


#Inputs
Testing = True #Set to false when used in ArcPro
Thread_Count = 1 # multiprocessing currently not setup correctly
# WSE_start = 511   # **********on the west Confluence of Erie canal/Genesse River**********
WSE_start = 511.6 # Eagle harbor
# WSE_start = 513 # for 60 mi 
WSE_slope = -1/316800
mile_marker_start = 294.6
# mile_marker_start = 320.6   # for 60 mi
# WSE_slope = -0.00000315394 # average normal pool slope
THRESH = 0.5 # threshold for checking crests are above invert and toe below invert


# In[3]:


# Fill in variables in StandaloneInputs function if testing
def StandaloneInputs():
    FilName = "pf.txt"
    # Output_Folder = "B:\\GIS_Data\\2021\\21P25022-070 - EAP Updates and Dam Breach Assessments\\Levee Script Testing\\Pilot Runs\\60 Mile"
    Output_Folder = "P:\\2022\\Albany\\21250022_070_Embankment Breach Modeling and Risk Assessment (Att. B)\\03-SEProducts\\07-GIS\\01-Data\\04-Shapefile\\Pilot Study\\Eagle Harbor Cross Section Data"
    Output_File_Name = "Results.csv"
    Output_File_Path = os.path.join(Output_Folder, Output_File_Name)
    Seg_ID = "NYPA"
    Final_Folder = Output_Folder
    # Final_Folder = "B:\\GIS_Data\\2021\\21P25022-070 - EAP Updates and Dam Breach Assessments\\Levee Script Testing\\Working_output"
    Section_spacing = 50
    return(Output_Folder+"\\"+FilName,Output_Folder,Output_File_Name,Output_File_Path,Seg_ID,Final_Folder,Section_spacing)

def ArcGISInputs():
    import arcpy
    FilName = arcpy.GetParameterAsText(0)
    Output_Folder = arcpy.GetParameterAsText(1)
    Output_File_Name = arcpy.GetParameterAsText(2)
    Output_File_Path = os.path.join(Output_Folder, Output_File_Name)
    Seg_ID = arcpy.GetParameterAsText(5)
    Final_Folder = arcpy.GetParameterAsText(4)
    Section_spacing = float(arcpy.GetParameterAsText(6))
    return(FilName,Output_Folder,Output_File_Name,Output_File_Path,Seg_ID,Final_Folder,Section_spacing) 


# In[4]:




Levee_StartTime = time.process_time() #Record the time at the start of script

# Gather Inputs
if Testing == False:
    FilName,Output_Folder,Output_File_Name,Output_File_Path,Seg_ID,Final_Folder,Section_spacing = ArcGISInputs()
else:
    FilName,Output_Folder,Output_File_Name,Output_File_Path,Seg_ID,Final_Folder,Section_spacing = StandaloneInputs()


# In[5]:



# Calculates relative maxes and mins
def min_max(num,df, n=40):
    n = 40  # number of points to be checked before and after
    dff = df[df['section']==num]
    dff['min'] = dff.iloc[argrelextrema(dff.z.values, np.less_equal,order=n)[0]]['z']
    dff['max'] = dff.iloc[argrelextrema(dff.z.values, np.greater_equal,order=n)[0]]['z']
    dfmin = dff[dff['min'].notnull()]
    dfmax = dff[dff['max'].notnull()]
    #dfmin = dfmin.drop(columns=['max','z'])
    #dfmax = dfmax.drop(columns=['min','z'])
    dfmin = dfmin.drop_duplicates(subset=['min'], keep='last')
    dfmax = dfmax.drop_duplicates(subset=['max'], keep='last')
    return dfmin, dfmax  


# In[6]:




######### Import and Organize Data

# Reads the text file as a Pandas dataframe
print("Reading data")
StartTime = time.process_time()
df = pd.read_fwf(FilName, header=None, names=['x','y','z'])
print("Time (min): ", (time.process_time()-StartTime)/60)

# Create path for output folders
XS_Folder = Output_Folder + "\\XS_Plots"
if not os.path.exists(XS_Folder):
    os.makedirs(XS_Folder)
Profile_Folder = Output_Folder + "\\Profile_Plots"
if not os.path.exists(Profile_Folder):
    os.makedirs(Profile_Folder)

# Seperate the imported file into individual x-sections
print("Setting individual sections")
StartTime = time.process_time()
indexes_to_drop=[]
df['section']=0
    
indexes_to_drop = df.index[df['x']=='END']    
indexes_to_keep = set(range(df.shape[0])) - set(indexes_to_drop)
df = df.take(list(indexes_to_keep))
indexes_Levee_change = df.index[df['z'].isnull()]

warnings.simplefilter(action='ignore')

for i in range(len(indexes_Levee_change)-1):
    df['section'].loc[indexes_Levee_change[i]:indexes_Levee_change[i+1]]=i+1

df['section'].loc[indexes_Levee_change[-1]:]=len(indexes_Levee_change)
df = df.dropna()

# set the end coordinates for each cross-section
df['x0'] = df.groupby('section')['x'].transform('first')
df['y0'] = df.groupby('section')['y'].transform('first')
df = df.astype({'x':float, 'y':float, 'z':float, 'section':int, 'x0':float, 'y0':float })
df['dx'] = df['x'] - df['x0']
df['dy'] = df['y'] - df['y0']
df['Xdata'] = np.sqrt(df['dx']**2 + df['dy']**2)

'''
unique_sectionbers = df['section'].unique()
df['x0'] = 0
df['y0'] = 0
        
for i in unique_sectionbers:
    x0 = df[df['section']==i].iloc[0]['x']
    y0 = df[df['section']==i].iloc[0]['y']
    df.x0.iloc[np.array(df['section']==i)] = x0
    df.y0.iloc[np.array(df['section']==i)] = y0
    
df = df.astype({'x':float, 'y':float, 'z':float, 'section':int, 'x0':float, 'y0':float }) 
print("Time (min): ", (time.process_time()-StartTime)/60)
# Calculating stationing in sections
print("Calculating stationing in sections")
StartTime = time.process_time()
df['dx'] = df['x']-df['x0']             
df['dy'] = df['y']-df['y0']            
# Calculate the cross-section 'Xdata' which is called stationing in the section
df['Xdata'] = np.sqrt(df['dx']**2 + df['dy']**2)
df.drop(['x', 'y', 'x0', 'y0', 'dx', 'dy'], axis=1, inplace=True)
df = df.reindex(columns=['x','y','z','section','Xdata','z'])
# df = df[(np.abs(stats.zscore(df['Xdata'])) < 3)]
print("Time (min): ", (time.process_time()-StartTime)/60)   
'''       
print(df)


# In[7]:



## Setup dataframe to hold results
columns_ = ['section','station','mile_post',
            'l_slop_out','l_slop_in','r_slop_in','r_slop_out',
            'l_toe_EL','r_toe_EL','l_crest_EL','r_crest_EL',
            'l_toe_station','r_toe_station','l_crest_station','r_crest_station',
           'invert', 'invert_station',
            'l_crest_w','r_crest_w','invert_w',
            'comp_time','link',
            'l_crest_E','l_crest_N',
            'r_crest_E','r_crest_N',
            'invert_E','invert_N',
            'l_toe_E','l_toe_N',
           'r_toe_E','r_toe_N', 
           'r_toe_hd', 'l_toe_hd',
           'r_height','l_height',
           'WSE', 'hyd_depth', 'avg_grad_R', 'avg_grad_L',
           'min_toe_L', 'min_toe_R', 'rep_slope_R', 'rep_slope_L','descending']
# results = pd.DataFrame(columns = columns_, index = [i for i in range(1,df['section'].nunique()+1)])
max_section = df['section'].nunique()
results = pd.DataFrame(columns=columns_, index=range(1, max_section + 1))


# In[53]:


## Function that processes one section at a time
def proc_sections(j):  
    n_refine=25   #min_max(j,df,n_refine)- if you write min_max(j,df) it consider the default value of 40 for n_refine
    ## This function should get moved 
    ##################################################HALFTRAPZ TURNED OFF
    def halfTrapz(xdata, per_1, per_2):
        # These variables should be added to the function call
        ydata = np.zeros(len(xdata))                      
        dfy = pd.DataFrame(list(zip(xdata,ydata)),columns=['x','y'])
        
        if (per_1 + per_2) > 1:
            per_sum = per_1 + per_2
            per_1 = per_1/per_sum
            per_2 = per_2/per_sum
        
        rng = xdata.iloc[-1] - xdata.iloc[0]
        m = (c-a)/(per_2*rng)
        
        #per_3 = 1 - per_1 - per_2
        
        point_count = len(xdata)
        
        ta1 = int(point_count*per_1)
        ta2 = int(point_count*(per_1 + per_2))

        # Calculate the y-coordinates
        dfy.iloc[:ta1,1] = a
        b = a - m*dfy.iloc[ta1,0] 
        dfy.iloc[ta1:ta2,1] =  m*dfy.iloc[ta1:ta2,0] + b 
        dfy.iloc[ta2:,1] = c
        return dfy['y'] 
        
    # Create section figure
    fig = plt.figure(0, figsize=(10,7))
    ax = fig.add_subplot(1, 1, 1)
       
    results['section'].iloc[j] = j
    results['station'].iloc[j] = j*Section_spacing
    results['mile_post'].iloc[j] = mile_marker_start + j*Section_spacing* (1/5280)   #converting ft to mi
    # Plot the water surface elevation
    WSE = results['station'].iloc[j]*WSE_slope + WSE_start
    ax.plot([0,df['Xdata'].max()], [WSE,WSE], '-b', label="WSE")
    
    print("Cross-section Number: ",j)
    StartTime = time.process_time() #Record time at start of current x-section

    # Calculate the relative maxes and mins
    dfmin, dfmax = min_max(j,df,n_refine) ############# why is this a function?

    # Plot the original data and the relative maxes and mins
    ax.plot(df[df['section']==j]['Xdata'], df[df['section']==j]['z'], 'k', label="Raw Data")
    ax.scatter(dfmin['Xdata'], dfmin['min'], s=40, c='r', marker='o', label="All Relative Mins")
    ax.scatter(dfmax['Xdata'], dfmax['max'], c='c', label="All Relative Maxs")

    # Calculate distance from the x-section center
    section_width = df['Xdata'].max()
    dfmin['cent_d'] = np.abs(dfmin['Xdata'] - section_width/2)
    
    success_bool = 1

    # Set the invert as the closest min to the center
    if len(dfmin)>0:
        invert = dfmin.iloc[dfmin['cent_d'].argmin()]
        invert_x = invert['Xdata']
        results['invert'].iloc[j] = invert['min']
        results['invert_station'].iloc[j] = dfmin.iloc[dfmin['cent_d'].argmin()]['Xdata']
        # Plot the invert
        ax.scatter(invert['Xdata'], invert['min'], s=100, c='k', marker='+', label="Invert")

        # Calculate the closest maxes to the invert
        dfmax['int_d'] = dfmax['Xdata'] - invert_x
        dfmax['int_d_abs'] = np.abs(dfmax['Xdata'] - invert_x)
        dfmax_right = dfmax[(dfmax['int_d']>0) & (dfmax['z'] > (invert['z'] + THRESH)) & (dfmax['Xdata'] < 0.95*section_width) & (dfmax['Xdata'] > 0.05*section_width)]
        dfmax_left = dfmax[(dfmax['int_d']<0)  & (dfmax['z'] > (invert['z'] + THRESH)) & (dfmax['Xdata'] < 0.95*section_width) & (dfmax['Xdata'] > 0.05*section_width)]
        
        results['invert_E'].iloc[j] = invert['x']
        results['invert_N'].iloc[j] = invert['y']
        results['WSE'].iloc[j] = WSE
        results['hyd_depth'].iloc[j] =  results['WSE'].iloc[j] - results['invert'].iloc[j]
           
        # if there is at least one max right of the invert
        if len(dfmax_right)>0:
            # find the closest max to the right of the invert
            dmax_right = dfmax_right.iloc[dfmax_right['int_d_abs'].argmin()]
            
            results['r_crest_E'].iloc[j] = dmax_right['x']
            results['r_crest_N'].iloc[j] = dmax_right['y']
            results['r_crest_EL'].iloc[j] = dmax_right['max']
            results['r_crest_station'].iloc[j] = dmax_right['Xdata']
            
            # plot the right crest
            ax.scatter(dmax_right['Xdata'], dmax_right['max'], s=100, c='m', marker='+', label="Northern Crest")
            
            # Curve fit invert to dmax_right
            nn = invert.name
            mm = dmax_right.name
            xdata = df[df['section']==j]['Xdata'].loc[nn:mm]
            ydata = df[df['section']==j]['z'].loc[nn:mm]
            a = invert['min']
            c = dmax_right['max']
            
            try:
                popt, pcov = curve_fit(halfTrapz, xdata, ydata, bounds=([0.0,0.0], [1.0, 1.0]), p0=[0.3,0.3], method='dogbox')

                per_1, per_2 = popt
                rng = xdata.iloc[-1] - xdata.iloc[0]
                slope = (c-a)/(rng*per_2)
                results['r_slop_in'].iloc[j] = np.abs(1/slope)
                Right_Crest_Left_Edge_xdata = xdata.iloc[0] + (per_1+per_2)*rng
                Right_Invert_xdata = xdata.iloc[0] + per_1*rng

                # Plot the best fit
                ax.plot(xdata, halfTrapz(xdata, *popt), 'r--')  

            except Exception:
                success_bool = 0
                print('inner right issue')
                traceback.print_exc()
                
            # Calculate the outboard toe
            right_toes = dfmin[(dfmin['Xdata']>dmax_right['Xdata']) & (dfmin['z'] < (dmax_right['z'] - THRESH)) & (dfmin['Xdata'] < 0.95*section_width) & (dfmin['Xdata'] > 0.05*section_width)]
                
            # if there is a min right of the right crest
            if len(right_toes)>0:
                right_toe = right_toes.iloc[right_toes['Xdata'].argmin()]

                # results['r_toe_E'].iloc[j] = right_toe['x']
                # results['r_toe_N'].iloc[j] = right_toe['y']
                # results['r_toe_EL'].iloc[j] = right_toe['min']
                # results['r_toe_station'].iloc[j] = right_toe['Xdata']

                # Plot the outboard toe
                # ax.scatter(right_toe['Xdata'], right_toe['min'], s=100, c='m', marker='x', label="Northern Toe")

                # Curve fit 
                nn = dmax_right.name
                mm = right_toe.name
                xdata = df[df['section']==j]['Xdata'].loc[nn:mm]
                ydata = df[df['section']==j]['z'].loc[nn:mm]
                a = dmax_right['max']
                c = right_toe['min']

                try:
                    popt, pcov = curve_fit(halfTrapz, xdata, ydata, bounds=([0.0,0.0], [1.0, 1.0]), p0=[0.3,0.3], method='dogbox')

                    per_1, per_2 = popt
                    rng = xdata.iloc[-1] - xdata.iloc[0]
                    slope = (c-a)/(rng*per_2)
 
                    results['r_slop_out'].iloc[j] = np.abs(1/slope)
                    Right_Crest_Right_Edge_xdata = xdata.iloc[0] + per_1*rng    
                    results['r_crest_w'].iloc[j] = Right_Crest_Right_Edge_xdata - Right_Crest_Left_Edge_xdata
                    Right_Toe_x = xdata.iloc[0] + (per_1 + per_2)*rng 
                    # results['r_toe_hd'].iloc[j] = WSE - results['r_toe_EL'].iloc[j]
                    # results['r_height'].iloc[j] = results['r_crest_EL'].iloc[j] - results['r_toe_EL'].iloc[j]

                    # Obtains index of the point closest to the calculated toe
                    index = int( (df[(df['section']==j)].Xdata - Right_Toe_x).abs().argsort()[:1])
                    
                    # Calculates the northing and easting at this index
                    # results['r_toe_E'].iloc[j] = df[(df['section']==j)].iloc[index].x
                    # results['r_toe_N'].iloc[j] = df[(df['section']==j)].iloc[index].y
                             
                    # Plot the best fit
                    ax.plot(xdata, halfTrapz(xdata, *popt), 'r--')

                except Exception:
                    success_bool = 0
                    print('inner right issue')
                    traceback.print_exc()
            else:
                success_bool = 0
        else:
            success_bool = 0

        
        # if there is at least one max left of the invert
        if len(dfmax_left)>0:         
            # find the closest max to the left of the invert
            dmax_left = dfmax_left.iloc[dfmax_left['int_d_abs'].argmin()]
            
            results['l_crest_E'].iloc[j] = dmax_left['x']
            results['l_crest_N'].iloc[j] = dmax_left['y']
            results['l_crest_EL'].iloc[j] = dmax_left['max']
            results['l_crest_station'].iloc[j] = dmax_left['Xdata']
            
            # plot the left crest
            ax.scatter(dmax_left['Xdata'], dmax_left['max'], s=100, c='g', marker='+', label="Southern Crest")

            # Curve fit dmax_left to invert
            nn = dmax_left.name
            mm = invert.name
            xdata = df[df['section']==j]['Xdata'].loc[nn:mm]
            ydata = df[df['section']==j]['z'].loc[nn:mm]
            c = invert['min']
            a = dmax_left['max']
   
            try:               
                popt, pcov = curve_fit(halfTrapz, xdata, ydata, bounds=([0.0,0.0], [1.0, 1.0]), p0=[0.3,0.3], method='dogbox')

                per_1, per_2 = popt
                rng = xdata.iloc[-1] - xdata.iloc[0]
                slope = (c-a)/(rng*per_2)
                results['l_slop_in'].iloc[j] = np.abs(slope)
                Left_Crest_Right_Edge_xdata = xdata.iloc[0] + (per_1)*rng 
                Left_Invert_xdata = xdata.iloc[0] + (per_1+per_2)*rng
                
                # Plot the best fit
                ax.plot(xdata, halfTrapz(xdata, *popt), 'r--')
                
            except Exception:
                success_bool = 0
                print('inner left issue')
                traceback.print_exc()
                
            # Calculate the outboard toe
            left_toes = dfmin[(dfmin['Xdata']<dmax_left['Xdata']) & (dfmin['z'] < (dmax_left['z'] - THRESH)) & (dfmin['Xdata'] < 0.95*section_width) & (dfmin['Xdata'] > 0.05*section_width)]               

            # if there is a min left of the left crest
            if len(left_toes)>0:
                left_toe = left_toes.iloc[left_toes['Xdata'].argmax()]
                # results['l_toe_E'].iloc[j] = left_toe['x']
                # results['l_toe_N'].iloc[j] = left_toe['y']
                # results['l_toe_EL'].iloc[j] = left_toe['min']
                # results['l_toe_station'].iloc[j] = left_toe['Xdata']

                # Plot the outboard toe
                # ax.scatter(left_toe['Xdata'], left_toe['min'], s=100, c='g', marker='x', label="Southern Toe")

                # Left toe to dfmax_left
                nn = left_toe.name
                mm = dmax_left.name
                xdata = df[df['section']==j]['Xdata'].loc[nn:mm]
                ydata = df[df['section']==j]['z'].loc[nn:mm]
                a = left_toe['min']
                c = dmax_left['max']

                try:                       
                    popt, pcov = curve_fit(halfTrapz, xdata, ydata, bounds=([0.0,0.0], [1.0, 1.0]), p0=[0.3,0.3], method='dogbox')

                    per_1, per_2 = popt
                    rng = xdata.iloc[-1] - xdata.iloc[0]
                    slope = (c-a)/(rng*per_2)
 
                    results['l_slop_out'].iloc[j] = np.abs(1/slope)
                    Left_Crest_Left_Edge_xdata = xdata.iloc[0] + (per_1 + per_2)*rng    
                    results['l_crest_w'].iloc[j] = Left_Crest_Right_Edge_xdata - Left_Crest_Left_Edge_xdata  
                    Left_Toe_x = xdata.iloc[0] + per_1*rng 
                    # results['l_toe_hd'].iloc[j] = WSE - results['l_toe_EL'].iloc[j]
                    # results['l_height'].iloc[j] = results['l_crest_EL'].iloc[j] - results['l_toe_EL'].iloc[j]

                    # Obtains index of the point closest to the calculated toe
                    # index = int( (df[(df['section']==j)].Xdata - Left_Toe_x).abs().argsort()[:1])
                    
                    # Calculates the northing and easting at this index
                    # results['l_toe_E'].iloc[j] = df[(df['section']==j)].iloc[index].x
                    # results['l_toe_N'].iloc[j] = df[(df['section']==j)].iloc[index].y
                    
                    # Plot the best fit
                    ax.plot(xdata, halfTrapz(xdata, *popt), 'r--')

                except Exception:
                    success_bool = 0
                    print('outter left issue')
                    traceback.print_exc()
            else:
                success_bool = 0
        else:
            success_bool = 0
                   
        if success_bool:
            results['invert_w'].iloc[j] = Right_Invert_xdata + Left_Invert_xdata        
        
    else:
        # no invert found
        success_bool = 0
        print("no data")
    
    # Set the tick marks and grid on plot
    major_ticks = np.linspace(0, 1000, num=11)
    minor_ticks = np.linspace(0, 1000, num=21)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.grid(which='both')
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    
    # Add labels and set axis limits to plot
    ax.set_xlim((0,df['Xdata'].max()))
    ax.set_xlabel('station (ft)')
    ax.set_ylabel('Elevation (ft)')
    ax.legend()
    Plot_title = f"Mile Post: {results['mile_post'].iloc[j]:.3f}"
    # Plot_title = str(Seg_ID) + " - " + str(j)
    plt.title(Plot_title)
    
    # Save the plot
    number = f"{results['mile_post'].iloc[j]:.3f}"
    new_miMarker = number.replace(".", "_")
    Save_name = f"Section_{j} MilePost {new_miMarker}.png"
    # Save_name = str(Seg_ID) + " - " + str(j) + '.png'
    results['link'].iloc[j] = XS_Folder + "\\" + Save_name
    
    # Record the processing time
    results['comp_time'].iloc[j] = (time.process_time()-StartTime)/60
    return (results, fig, ax)


# In[97]:


'''
Definitions: *Please use these terms to be more readable*
cond_I: condition for Invert
cond_RC: condition for Right Crest
cond_RT: condition for Right Toe
cond_LC: condition for Left Crest
cond_LT: condition for Left Toe
cond_WSE: condition for water surface elevation
inner_L: Left side wall 
inner_R: Right side wall
'''
def gradient_calc(df, results,j, fig, ax, which_grad='Right'):
# def gradient_calc(df, results,j, fig, ax, left_grad=True, right_grad=True):
    
    # gradient = pd.DataFrame()
    cond_I = results.invert_station.iloc[j] 
    cond_WSE = results.WSE.iloc[j]
    WSE = results.WSE[results.section==j].values
    max_grad = 0  
    y1=[];y2=[];x1=[];x2=[] 
    print(which_grad)
    
    if which_grad=='Right':  
        cond_RC = results.r_crest_station.iloc[j]
        inner_R = df[(df.section==j) & (df.Xdata>= cond_I) & (df.Xdata< cond_RC) & (df.z < cond_WSE)]['Xdata'].values
        inner_general = inner_R
        y2_val = df[(df.section==j) & (df.Xdata> cond_RC) & (df.z < cond_WSE)]['z'].values # n array of outer y values
        x2_val = df[(df.section==j) & (df.Xdata> cond_RC) & (df.z < cond_WSE)]['Xdata'].values # n array of outer x values         
    elif which_grad=='Left':
        cond_LC = results.l_crest_station.iloc[j]
        inner_L = df[(df.section==j) & (df.Xdata<= cond_I) & (df.Xdata> cond_LC) & (df.z < cond_WSE)]['Xdata'].values
        inner_general = inner_L
        y2_val = df[(df.section==j) & (df.Xdata< cond_LC) & (df.z < cond_WSE)]['z'].values # n array of outer y values
        x2_val = df[(df.section==j) & (df.Xdata< cond_LC) & (df.z < cond_WSE)]['Xdata'].values # n array of outer x values         
    else:
        print('PLEASE TYPE EITHER "Right" OR "Left" WITH EXACT CAPs')
        # break

    
    for x1_val in inner_general:
        y1_val = df[(df.section == j) & (df.Xdata == x1_val)]['z'].item() # scalar of y value

        length_grad_val = (np.sqrt((y2_val-y1_val)**2 + (x2_val-x1_val)**2)) # n array of distances
        del_head = WSE - y2_val
        
        gradients = (del_head/length_grad_val) * (np.abs(x2_val-x1_val)/length_grad_val)


        if len(gradients)==0:
            return (df, results,j, fig, ax,which_grad)        

        gradient = max(gradients)
        if gradient > max_grad:
            max_grad = gradient
            i = gradients.argmax()
            x1 = x1_val
            y1 = y1_val
            y2 = y2_val[i]
            x2 = x2_val[i]
    

    if (np.isnan(x1) ==False) & (np.isnan(y1)==False) & (np.isnan(x2)==False) & (np.isnan(y2)==False):
        ax.plot([x1, x2], [y1, y2], label=f'Max Gradient_{str(which_grad)[0]}: {gradient:.2f}')
        if which_grad =='Right':
            results['avg_grad_R'].iloc[j] = gradient
            # results['min_toe_R'].iloc[j] = min(y2, results.r_toe_EL.iloc[j]) 
            ######
            WSE_EL_R = df[(df.section == j) & (df.z<cond_WSE) & (df.Xdata> cond_RC) & (df.Xdata< x2)]['z']
            if not WSE_EL_R.empty:
                WSE_EL_R = WSE_EL_R.sort_values(ascending=False).iloc[0] 
                WSE_station_R = df[(df.section == j) & (df.z==WSE_EL_R)]['Xdata'].iloc[0]
                firstX = WSE_station_R
                lastX = x2
                firstY = WSE_EL_R
                lastY=y2
                slope_WSE_gradient_R = np.abs((lastY-firstY)/(lastX-firstX))
                results['rep_slope_R'].iloc[j] = 1/slope_WSE_gradient_R
                ax.plot([firstX,lastX],[firstY, lastY], label=f'rep_slope_N: {1/slope_WSE_gradient_R:.2f}H:1V', ls='-.', color='gold') 
                
                results['r_toe_station'].iloc[j] = lastX
                results['r_toe_EL'].iloc[j] = lastY
                
                results['r_toe_E'].iloc[j] = df[(df.Xdata == lastX) & (df.section==j)]['x'].values[0]
                results['r_toe_N'].iloc[j] = df[(df.z == lastY) & (df.section==j)]['y'].values[0]
                # left_toe = left_toes.iloc[left_toes['Xdata'].argmax()]
                # results['l_toe_E'].iloc[j] = left_toe['x']
                # results['l_toe_N'].iloc[j] = left_toe['y']
                ax.scatter(lastX, lastY, s=100, c='m', marker='x', label="Northern Toe")
                
                results.r_height.iloc[j] = results.r_crest_EL.iloc[j] - results.r_toe_EL.iloc[j]
                if results.r_height.iloc[j] < 200:
                    results.r_toe_hd.iloc[j] = results.WSE.iloc[j] - results.r_toe_EL.iloc[j]
            ########################################################
                    

        elif which_grad =='Left':
            results['avg_grad_L'].iloc[j] = gradient
            results['min_toe_L'].iloc[j] = min(y2, results.l_toe_EL.iloc[j])
            
            WSE_EL_L = df[(df.section == j) & (df.z<cond_WSE) & (df.Xdata< cond_LC) & (df.Xdata> x2)]['z']
            if not WSE_EL_L.empty:  
                WSE_EL_L = WSE_EL_L.sort_values(ascending=False).iloc[0]
                WSE_station_L = df[(df.section == j) & (df.z==WSE_EL_L)]['Xdata'].iloc[0]
                firstX = x2
                lastX = WSE_station_L
                firstY = y2
                lastY = WSE_EL_L
                slope_WSE_gradient_L = np.abs((lastY-firstY)/(lastX-firstX))
                results['rep_slope_L'].iloc[j] = 1/slope_WSE_gradient_L
                ax.plot([firstX,lastX],[firstY, lastY], label=f'rep_slope_S: {1/slope_WSE_gradient_L:.2f}H:1V', ls='-.', color='pink') 
    
                results['l_toe_station'].iloc[j] = firstX
                results['l_toe_EL'].iloc[j] = firstY
                
                results['l_toe_E'].iloc[j] = df[(df.Xdata == firstX) & (df.section==j)]['x'].values[0]
                results['l_toe_N'].iloc[j] = df[(df.z == firstY) & (df.section==j)]['y'].values[0]
                ax.scatter(firstX, firstY, s=100, c='b', marker='x', label="Southern Toe")
                
                results.l_height.iloc[j] = results.l_crest_EL.iloc[j] - results.l_toe_EL.iloc[j]
                if results.l_height.iloc[j] < 200:
                    results.l_toe_hd.iloc[j] = results.WSE.iloc[j] - results.l_toe_EL.iloc[j]
    # cond_RC = results.r_crest_station.iloc[j]
    # cond_LC = results.l_crest_station.iloc[j] 
    # left_desc = df[(df.section == j) & (df.Xdata < cond_LC)]['z'].reset_index()
    # right_desc = df[(df.section == j) & (df.Xdata > cond_RC)]['z'].reset_index()
    
    # fill_in_WSE = df[(df.section == j) & (df.z<=cond_WSE) & (df.Xdata<cond_RC) & (df.Xdata>cond_LC)]   
    # ax.fill_between(fill_in_WSE['Xdata'], fill_in_WSE['z'], cond_WSE, where=(fill_in_WSE['z'] < cond_WSE), color='blue', alpha=0.5)
    return (df, results,j, fig, ax,which_grad)




def descending(df, results, j, side, diff=2):
    col ='z'
    cond_RC = results.r_crest_station.iloc[j]
    cond_LC = results.l_crest_station.iloc[j] 
    left_desc = df[(df.section == j) & (df.Xdata < cond_LC)]['z'].reset_index()
    right_desc = df[(df.section == j) & (df.Xdata > cond_RC)]['z'].reset_index()
    
    if side.upper() == 'LEFT':
        operator = '>'
        dfn=left_desc 
    elif side.upper() == 'RIGHT':
        operator = '<'
        dfn=right_desc
        
    descending = True
    for i in range(len(dfn)):
        for j in range(i, len(dfn) - 1):
            if operator == '>' and dfn[col].iloc[i] > dfn[col].iloc[j + 1]:
                break
            elif operator == '<' and dfn[col].iloc[i] < dfn[col].iloc[j + 1]:
                break
            else:
                if np.abs(dfn[col].iloc[j + 1] - dfn[col].iloc[i]) < diff:
                    continue
                else:
                    descending = False
                    # results.iloc[j]=descending
                    return descending
    
    # results.iloc[j]=descending                
    return descending


# In[99]:


def figs(fig, ax):
            
    '''
    # For 1:1 plotting
    plt.axis('equal')
    plt.xlim([475,650])
    plt.ylim([497,515])
    '''
    
    ax.legend()
    
    # Create path for output folders
    XS_Folder = Output_Folder + "\\XS_Plots"
    if not os.path.exists(XS_Folder):
        os.makedirs(XS_Folder)
    Profile_Folder = Output_Folder + "\\Profile_Plots"
    if not os.path.exists(Profile_Folder):
        os.makedirs(Profile_Folder)
    # Save the plot
    
    
    number = f"{results['mile_post'].iloc[j]:.3f}"
    new_miMarker = number.replace(".", "_")
    Save_name = f"Section_{j} MilePost {new_miMarker}.png"
    # Save_name = str(Seg_ID) + " - " + str(j) + '.png'
    results['link'].iloc[j] = XS_Folder + "\\" + Save_name
    fig.savefig(XS_Folder + "\\" + Save_name,dpi=300)    
    
    # plt.figure().clear()
    # plt.close()
    # plt.cla()
    # plt.clf() 
    return 


# In[12]:

############################### Just a test to know the cross-sections of natuaral drainage areas.
# bad_cross=[]
# for j in range(1,df['section'].nunique()):
#         # results, fig, ax = proc_sections(j)
#         # plt.figure().clear()
#         # plt.close()
#         # plt.cla()
#         # plt.clf()
#         right_desc = descending(df, results, j, 'Right', diff=0.5)
#         left_desc = descending(df, results, j, 'Left', diff=0.5)

#         if (right_desc==True) | (left_desc==True):
#             results['descending'].iloc[j]='X'
#             print(f'check sec {j}')
#             bad_cross.append(j)
#             # results, fig, ax = proc_sections(j)
#             # df, results,j, fig, ax, which_grad = gradient_calc(df, results,j, fig, ax, 'Right')
#             # # right_desc = descending(df, results, j, 'Right', diff=0.5)
#             # figs(fig, ax)
#             # df, results,j, fig, ax, which_grad = gradient_calc(df, results,j, fig, ax, 'Left')
#             # figs(fig, ax)
#             # # left_desc = descending(df, results, j, 'Left', diff=0.5)
#             # plt.figure().clear()
#             # plt.close()
#             # plt.cla()
#             # plt.clf()




## Process the Sections
if Thread_Count ==1:
    ## Run the cross-sectional analysis in series
    
    # for j in range(4693,4694):
    for j in range(1,df['section'].nunique()):
        results, fig, ax = proc_sections(j)
        df, results,j, fig, ax, which_grad = gradient_calc(df, results,j, fig, ax, 'Right')
        # right_desc = descending(df, results, j, 'Right', diff=0.5)
        figs(fig, ax)
        df, results,j, fig, ax, which_grad = gradient_calc(df, results,j, fig, ax, 'Left')
        figs(fig, ax)
        # left_desc = descending(df, results, j, 'Left', diff=0.5)
        plt.figure().clear()
        plt.close()
        plt.cla()
        plt.clf()
        # if (right_desc==True) | (left_desc==True):
        #     results['descending'].iloc[j]=='X'
        #     print(f'check sec {j}')
        
 

        '''            
            plt.figure().clear()
            plt.close()
            plt.cla()
            plt.clf()
            pass
        '''
else:
    ## Run the cross-sectional analysis in parallel    
    with mp.Pool(Thread_Count) as p:
        p.map(proc_sections, [i for i in range(1,df['section'].nunique())])
 


# In[105]:



### Save Summary Data   
Summary_Data_File_Path = os.path.join(Final_Folder, "Summary_Results.csv")
results.to_csv(Summary_Data_File_Path)

# Export results files for non null data
results[results['invert_E'].notnull()].to_csv(os.path.join(Final_Folder, "Invert.csv"))
results[results['l_crest_E'].notnull()].to_csv(os.path.join(Final_Folder, "Southern_Crest.csv"))
results[results['r_crest_E'].notnull()].to_csv(os.path.join(Final_Folder, "Northen_Crest.csv"))
results[results['l_toe_E'].notnull()].to_csv(os.path.join(Final_Folder, "Southern_Toe.csv"))
results[results['r_toe_E'].notnull()].to_csv(os.path.join(Final_Folder, "Northen_Toe.csv"))
# bad_cross.to_csv
#### Summary of a summary table that Allan asked:
# columns_to_report = results[['section','station','mile_post','l_crest_w','r_crest_w','r_height','l_height', 'invert','hyd_depth', 'WSE','rep_slope_R', 'rep_slope_L']]
# columns_to_report.dropna(axis=0, how='all', inplace=True)            
# columns_to_report.set_index(['section'], inplace=True) 
# # columns_to_report['hyd_head']=columns_to_report['WSE']-columns_to_report['invert']
# columns_to_report.to_csv(os.path.join(Final_Folder, "Summary_Report_Allan.csv"))


### Plot Levee Profile with outliers (Plot for checking)   
plt.plot(results['mile_post'], results['invert'], 'k', label="Invert")
plt.plot(results['mile_post'], results['l_toe_EL'], '--r', label="Southern toe")
plt.plot(results['mile_post'], results['r_toe_EL'], '--g', label="Northern toe")
plt.plot(results['mile_post'], results['l_crest_EL'], '-r', label="Southern Crest")
plt.plot(results['mile_post'], results['r_crest_EL'], '-g', label="Northern Crest") 
plt.plot(results['mile_post'],  results['WSE'], '--b', label="WSE") 

lgd = plt.legend(loc="upper left", bbox_to_anchor=(1,1), borderaxespad=1)
plt.xlabel('mile_post')
plt.ylabel('Elevation (ft)')

Plot_title = str(Seg_ID)
plt.title(Plot_title) # Creates a plot title
Save_name = Profile_Folder + "\\" + str(Seg_ID) + ' - Profile With Outliers.png'
plt.savefig(Save_name,dpi=400,bbox_extra_artists=(lgd,), bbox_inches='tight') # Saves the image

plt.clf() # Clears the plot  
plt.close()

# results['WSE'] = -np.array(results['station'])*WSE_slope + WSE_start



fig = px.line(results, x='mile_post', y=['min_toe_L', 'l_crest_EL', 'invert', 'l_crest_EL', 'r_crest_EL','min_toe_R', 'WSE'])
px.line()
# fig.write_html(Profile_Folder + "\interactive_profile.html")


html = fig.to_html(full_html=False)
html = html.encode("UTF-8")
with open(Profile_Folder + "\interactive_EageleHarbor_profile.html", "wb") as f:
    f.write(html)
### Plot Levee Profile without Outliers (Clean plot)


# In[ ]:



### Output the file into the ArcPro model
if Testing==False: 
    import arcpy
    arcpy.SetParameter(3, Output_File_Path)


# In[ ]:




