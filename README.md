# Levee-optimization-CrossSection

This Python script performs cross-sectional analysis on a topographic cross-sections cut at a standard spacing along a canal centerline. The script calculates several geometric properties, such as crest elevations, channel invert, and toe elevations.

This version is for **Eagle Harbor**. To change it for other locations please update the following parameters:
>**a)** Water Surface Level: WSL need to be updated for eagle horbor we put 511.6 on the left (west) side and it decrease with the slop of 1/316800  
**b)** Mile posts starts at 296.5 on the left (west) side. For other locations it need to be updated.

### These are the following parameters:
> **section:** Section numbers.  
**station:**  
**mile marker:** mile marker.              
**l_slop_out:** Left side outboard.  
**l_slop_in:** Left side inboard.  
**r_slop_in:** Right side inboard.  
**r_slop_out:** Right side outboard.  
**l_toe_EL:** Left toe elevation.  
**r_toe_EL:** Right toe elevation.  
**l_crest_EL:** Left crest elevation.      
**r_crest_EL:** Right crest elevation.      
**l_toe_station:** Left toe station.  
**r_toe_station:** Right toe station.  
**l_crest_station:** Left crest station.  
**r_crest_station:** Right crest station.  
**invert:** Invert.  
**invert_station:** Right station.  
**l_crest_w:** Left crest width.  
**r_crest_w:** Right crest width.  
**invert_w:** Invert width.  
**r_height:** WSL to right toe.  
**l_height:** WSL to left toe.  
**WSE:** Water Surface Elevation.  
**hyd_depth:** Hydraulic depth.  
**avg_grad_R:** Max gradient on right slope.  
**avg_grad_L:** Max gradient on left slope.  
**min_toe_L:** Minimum elevation between toe and gradient outlet elevation on left side.  
**min_toe_R:** Minimum elevation between toe and gradient outlet elevation on right side.  
**Rep. slope_R:** Slope from the WSE to gradient outlet on right side.   
**Rep. slope_L:** Slope from the WSE to gradient outlet on left side.  
