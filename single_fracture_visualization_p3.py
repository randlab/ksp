#C:\Users\A480325\miniconda\python.exe -m pip install statistics
#RUNS WITH PYHON 2.7
# C:\Users\A480325\miniconda\python.exe 
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt  #plt plt plt
from matplotlib import cm
from matplotlib.patches import Ellipse
import matplotlib.colors as colors
import math
import statistics
from matplotlib.collections import EllipseCollection
from PIL import Image
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from PIL import Image
from PIL import ImageDraw

from tqdm import tqdm  #the progress bar

from datetime import datetime
import termcolor
import os
#os.system('color')
import time

from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
import matplotlib.patches as mpatches

#print (termcolor.colored ('\nSTART %s\n','magenta')) % datetime.now()
print ('\nSTART %s' % datetime.now())

print (termcolor.colored ('\nreading coordinates and list of years with data...','cyan')) 

input_filepath = 'C:/Daten/KARST_Nagra/Dreybrodt/percolation_network/export/1K.txt'
#input_filepath = 'C:/Daten/KARST_Nagra/Dreybrodt/percolation_network/export/high_dip_04_WITH_initial_report.txt' #CASE 1 with initial report version
#input_filepath = 'C:/Daten/KARST_Nagra/Dreybrodt/percolation_network/export/high_dip_05_WITH_initial.txt' # CASE 2 with initial report version
#input_filepath = 'C:/Daten/KARST_Nagra/Dreybrodt/percolation_network/export/high_dip_07_refined.txt' # CASE 3 
#input_filepath = 'C:/Daten/KARST_Nagra/Dreybrodt/percolation_network/export/APR_2021/highdip_09_6_nodes_C.txt' # APR_2021
#input_filepath = 'C:/Daten/KARST_Nagra/Dreybrodt/percolation_network/export/APR_2021/case04_125K_50m_B.txt' # APR_2021
#input_filepath = 'C:/Daten/KARST_Nagra/philippe/fracture_networks/Case_4/high_dip_09_laminar.txt'#input_filepath = 'Y:/Feflow/export/1K.txt'
#nput_filepath = 'C:/Daten/KARST_Nagra/Dreybrodt/percolation_network/export/highdip_09_11nodes_50m.txt'
#input_filepath = 'C:/Daten/KARST_Nagra/philippe/fracture_networks/export/1000x2000_03_debug_version.txt'
#input_filepath = 'C:/Daten/KARST_Nagra/philippe/fracture_networks/export/case_02_100m_1e-3_initial.txt'   # _initial


with open(input_filepath) as f:
    content = f.readlines()
    
#look for coordinates block
data_export_idx = []

for line in content:
    if 'fracture' in str(line):
        coor_start = content.index(line) + 1 
        print ('\ncoor_start line : %s' % coor_start)
    if 'Data export' in str(line):
        data_export_idx.append(content.index(line)+1)
        
coor_end = data_export_idx[0] - 2
print ('coor_end line : %s' % coor_end)

coor_list= content[coor_start:coor_end]
coor_number = coor_end-coor_start
print ('coor_number : %s'% coor_number)

#store coordinates in array
coor_array = np.zeros ((coor_number,5)) # Y, X

for index,line in enumerate(coor_list):
    array_row = line.split(',')
    #print (len(array_row))
    for idx,value in enumerate(array_row):
        #print ('value : %s' % value)
        coor_array[index,idx] = value
        
upper_y = 500
lower_y = 400
filtered_coords_counter = 0

for index,line in enumerate(coor_array):
    if ( (line[2] <= upper_y) and (line[4]>=lower_y) ):
        filtered_coords_counter = filtered_coords_counter + 1

filtered_coords = np.zeros ((filtered_coords_counter,5))
filtered_coords_counter = 0

for index,line in enumerate(coor_array):
    if ( (line[2] <= upper_y) and (line[4]>=lower_y) ):
        for idx,value in enumerate(line):
            filtered_coords [filtered_coords_counter,idx] = value
        filtered_coords_counter = filtered_coords_counter + 1
    
# get list of years with data_export_idx
years =[]

for idx in data_export_idx:
    line_temp = content[idx]
    newstr = ''.join((ch if ch in '0123456789.' else ' ') for ch in line_temp)
    years.append ( [float(i) for i in newstr.split()] )
    
outflows = np.zeros ((len(years),1))
outflow_fractures = [999]   #case_01 1352 case02 = 5764 high_dip_09 = [6578,6579,6753] 38814,38897     [17093,17094,17138,17139,17140,17178,17179,17180]
outlet_coords = [1000,0.5] #case_01 50.4098 case 2 145  2000,50

    
for loop_idx,timestep in enumerate(years):
    #print ('loop_idx : %s'  % loop_idx)
    year = 'time = ' + str(timestep[0]) 
    #year = 'time = 10020 years'
    #print (year)

    print ( (termcolor.colored ('\n(%s/%s) plot made with %s years','white')) % (loop_idx+1,len(years),year) )

    for line in tqdm(content):
        if year in str(line):
            data_start = content.index(line) + 2
            data_end = data_start + coor_number
    
    if 'data_start' in locals():   
        #print ('')
        pass
    else:
        print (termcolor.colored ('\nyear not found','red')) 
        exit()        
            
    #print ('data_start : %s' % data_start)
    #print ('data_end : %s ' % data_end)
    
    #data_array
    
    data_list= content[data_start:data_end]
    data_array = np.zeros ((coor_number,6)) # Y, X
    
    for index,line in enumerate(data_list):
        array_row = line.split(',')
        #print (len(array_row) )
        for idx,value in enumerate(array_row):
            #print ('value : %s' % value)
            data_array[index,idx] = value   
    
    #PLOT      
    
    #colorbar normalization
    color_data_min = data_array[:,1].min()  # min value for colorbar
    color_data_max = data_array[:,1].max()  # max value for colorbar
    print ('Aperture / color data min : %s / color data max : %s' % (color_data_min, color_data_max) )
    
    #norm = mpl.colors.Normalize(vmin=color_data_min, vmax= color_data_max)
    norm = mpl.colors.LogNorm(vmin=color_data_min, vmax= color_data_max)
    #norm = mpl.colors.PowerNorm(vmin=color_data_min, vmax= color_data_max,gamma=2)
    
    #for colorbar normalization
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet) #JET
    #cmap.set_array([])
     
    print ( (termcolor.colored ('(%s/%s)lines collection generation and plotting...','cyan')) % (loop_idx+1,len(years)) )
    
    segs = []
    colors = []
    heads = np.zeros ((coor_number+1,3))
    color_flow_eq = []
    
     
    for index,row in enumerate(tqdm(coor_array)):
        #print ('index %s'% index)
        #print (row)
        x1 = coor_array[index,1]
        x2 = coor_array[index,3]
        y1 = coor_array[index,2]
        y2 = coor_array[index,4]
        
#        heads[index,0] = (x2+x1)/2
#        heads[index,1] = (y2+y1)/2 
#        heads[index,2] = data_array[index,6]
        
        color_flow_eq.append(cm.viridis(data_array[index,2]-1))
        
        #outflow plot
        cum_outfrac_flow = 0  #outflow plot
        for out_frac in outflow_fractures:
            if index == out_frac:
                cum_outfrac_flow = cum_outfrac_flow + data_array[index,4]
               # print (cum_outfrac_flow)
        
        #color_value = (data_array[index,1]-color_data_min)/(color_data_max-color_data_min) 
    
        color_value = data_array[index,1] * (1 / color_data_min) # log
        color_value = np.log10(color_value) /( np.log10(color_data_max)-np.log10(color_data_min)) # log
        
        #alpha = 2
        #color_value =  pow(data_array[index,1],alpha) /(pow(color_data_max,alpha)-pow(color_data_min,alpha))    
            
        #print (data_array[index,1])
        #print ('color value : %s' % color_value)   
        colors.append(cm.jet(color_value))
        segs.append( ( (x1, y1) , (x2, y2) ) )
    
    outflows[loop_idx] = cum_outfrac_flow  #outflow plot
    
    print ( (termcolor.colored ('(%s/%s)Fracture plot generation ( image 1 )','cyan')) % (loop_idx+1,len(years)) )
    

    #ln_coll = LineCollection(segs, colors=colors,linewidth = 1.5)
    ln_coll = LineCollection(segs, colors=colors,linewidth = 20)
   
    ax = plt.gca()
    
    #heads[-1,0] = 0
    #heads[-1,1] = 0
    #heads[-1,2] = 120 #case 1 50, case 2 145
  
# No heads in old files    
#    cont = ax.tricontour (heads[:,0],heads[:,1],heads[:,2],10, 
#        colors='grey', linestyles='dashed',linewidths=0.75)
#    ax.clabel (cont, fmt='%1.1f', fontsize=8,colors='black', inline=True)
    
    ax.add_collection(ln_coll)
    ax.set_xlim(0, max(coor_array[:,1])+ 10 )    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    #ax.set_ylim(0.4, 0.6)
    #ax.set_ylim(0.9* min(coor_array[:,2]), max(coor_array[:,2]) )
    #ax.set_ylim(0.9* min(coor_array[:,2]), 500+10 )
    ax.set_ylim(0, 1 )
    real_year = timestep[0] 
    real_year = '{:.1f}'.format(real_year)
    print ('Plot title / real year : %s' % real_year )
    plt.title('T = '+ str(real_year) + ' years')
    plt.xlabel('X [m]',labelpad = 10)
    plt.ylabel('Y [m]',labelpad = 10)
    
        #OUTLETS AND INLETS in fracture width plot
   
    
    #ax.plot(outlet_coords[0],outlet_coords[1], '>',color='red',markersize=5, label='outlets') #OUTLET
    ax.plot(outlet_coords[0]+1,outlet_coords[1], '>',color='red',markersize=5, label='outlets') #OUTLET
    #for idx,node_y in enumerate(coor_array[:,2]):  
        #if node_y == 500:
        #    plt.plot(coor_array[idx,1],node_y, 'v', color='limegreen', markersize=5, label='Inlets') #Inlets
        
#    nodes = [0,200,400,600,800,1000,1200,1400,1600,1800]
#    for node in nodes:
#        plt.plot(node,500, 'v', color='limegreen', markersize=5, label='Inlets') #Inlets
#        plt.plot(node,0, 'v', color='limegreen', markersize=5, label='Inlets') #Inlets
    
    plt.plot(-1,0.5, '>', color='limegreen', markersize=5, label='Inlets') #Inlets    
    
    outlets_legend = Line2D([0], [0],marker='>', color='red', markersize=5, label='Outlets',linewidth=0)
    inlets_legend = Line2D([0], [0],marker='>', color='limegreen', markersize=5, label='Inlets',linewidth=0)
            
    legend_1 = plt.legend(handles=[outlets_legend,inlets_legend],loc='best', 
        ncol = 2, bbox_to_anchor=(-0.2, -0.55, 0.5, 0.5) )
    
    #legend ISOLINES
#    isolines_legend = Line2D([0], [0], color='grey', linestyle = 'dashed', 
#        linewidth = 0.75, label ='Hydraulic head isolines [m]')
#    legend_2  = plt.legend(handles=[isolines_legend],loc='best', ncol = 2,
#        bbox_to_anchor=(0.35, -0.55, 0.5, 0.5) )
#        
#    plt.gca().add_artist(legend_1)
#    plt.gca().add_artist(legend_2)
    

    #colorbar    
    clb = plt.colorbar(cmap)#,format='%.0e')#, ticks=clb_ticks )
    clb.set_label('Aperture [m]',labelpad=-40, y=1.05, rotation=0)
        
    #for initial plot    
    clb.ax.tick_params(labelsize=9)
    clb.ax.set_yticklabels(["{:.0e}".format(i) for i in clb.get_ticks()])
    
    figure = plt.gcf()
    figure.set_size_inches(20, 6)
    
    export_file = str(loop_idx) + '_' + str(real_year) + '_fracs.png'
    plt.savefig(export_file,dpi=300,bbox_inches = "tight")
    plt.clf()
    plt.close('all')
    
    #-------------------------------------------------------------------
    print ( (termcolor.colored ('(%s/%s)Flow Equation plot generation ( image 3 )','cyan')) % (loop_idx+1,len(years)) )
    #LAMINAR OR TURBULENT
    ax = plt.gca()
    
    line_flow_eq = LineCollection(segs, colors=color_flow_eq,linewidth = 20)
#   line_flow_eq = LineCollection(segs, colors=color_flow_eq,linewidth = 1.5
    ax.add_collection(line_flow_eq)
    ax.set_xlim(0, max(coor_array[:,1])+ 10 )    
    ax.plot(outlet_coords[0],outlet_coords[1], "or") #OUTLET
    #ax.set_ylim(0.9* min(coor_array[:,2]), 500+10 )
    ax.set_ylim(0, 1 )
    real_year = timestep[0] 
    real_year = '{:.1f}'.format(real_year)
    plt.title('T = '+ str(real_year) + ' years')
    plt.xlabel('X [m]',labelpad = 10)
    plt.ylabel('Y [m]',labelpad = 10)
    
    patch_1 = mpatches.Patch(color= cm.viridis(0), label='Laminar Flow')
    patch_2 = mpatches.Patch(color= cm.viridis(0.9999), label='Turbulent Flow')
    plt.legend(handles=[patch_1,patch_2],loc='best', ncol = 2,
        bbox_to_anchor=(0.35, -0.55, 0.5, 0.5) )
    
    figure = plt.gcf()
    figure.set_size_inches(20, 6)
    
    export_file = str(loop_idx) + '_' + str(real_year) + '_flow_eq.png'
    plt.savefig(export_file,dpi=300,bbox_inches = "tight")
    plt.clf()
    plt.close('all')
    
    #-------------------------------------------------------------------
    #Flow_evolution
    print ( (termcolor.colored ('(%s/%s)Flow Rate plot generation ( image 4 )','cyan')) % (loop_idx+1,len(years)) )
    
    outflows_liters = outflows / 86.4
    
    print (data_start)
    print (data_end)
    print (cum_outfrac_flow)
    #print (loop_idx)
    #number = outflows_liters[loop_idx]
    #print (np.format_float_scientific(number,precision = 2, exp_digits=3))
    
    years_for_flow = years[0:len(outflows)]
    plt.scatter(years_for_flow,outflows_liters,s=5)  #m3/day to L/s
    #for line in outflows:
    #    print(line)
    
    max_years = str(max(years))
    max_years = float(max_years[1:-1])

    plt.xlim( ( 0, max_years )) 
    #plt.xlim( ( 0, 1481.71 )) 
    plt.ylim((1e-4,6000))
    #plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Simulation years', fontsize=10,labelpad=10)
    plt.ylabel('Outflow [L/s]', fontsize=10,labelpad=10)    
    plt.title('T = '+ str(real_year) + ' years')
    plt.grid(b=True, axis='x', which='major', color='black', linestyle='-',linewidth=0.2,alpha=0.5)
    plt.grid(b=True, axis='y', which='major', color='black', linestyle='-',linewidth=0.2,alpha=0.5)    
    
    figure = plt.gcf()
    figure.set_size_inches(4, 4)
    export_file = str(loop_idx) + '_' + str(real_year) + '_outflow.png'
    plt.savefig(export_file,dpi=150,bbox_inches = "tight")
    
    plt.clf()
    plt.close('all')
    
    #**********************************************************************
    #HISTOGRAM aperture
    
    print ( (termcolor.colored ('(%s/%s)Fracture histogram generation ( image 2 )','cyan')) % (loop_idx+1,len(years)) )
    
    log_bins = np.logspace(np.log10(5e-5),np.log10(5),100)
    
    filtered_apertures = np.zeros ((filtered_coords_counter,1))
    filtered_coords_counter = 0
    
    for index,line in enumerate(data_array):
        if line[0] in filtered_coords[:,0]:
            filtered_apertures[filtered_coords_counter] = line[1]
            filtered_coords_counter = filtered_coords_counter + 1
    
    #print ( len(filtered_apertures) )   
       
    #vertical seletion of fractures!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    #n,bins, patches = plt.hist(filtered_apertures,bins=log_bins)
    n,bins, patches = plt.hist(data_array[:,1],bins=log_bins)
    idx = 0
    
    for c, p in zip(bins, patches):    
        color_value = log_bins[idx] * (1 / color_data_min) # log
        color_value = np.log10(color_value) /( np.log10(color_data_max)-np.log10(color_data_min)) # log
        #print ('idx: %s color value : %f' % (idx, color_value))
        #time.sleep(0.5)
        plt.setp(p, 'facecolor', cm.jet(color_value))
        idx+=1
    
    #plt.hist(data_array[:,1],density=True, bins=200)
    plt.xlim((0.00005,50)) #plt.xlim((0.00005,50))
    plt.ylim((0.5,15000))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Aperture [m]', fontsize=10,labelpad=10)
    plt.ylabel('Frequency', fontsize=10,labelpad=10)
    
    ax = plt.gca()
    y_labels = ax.get_yticks()
    ax.set_yticklabels(['{:,.0f}'.format(x) for x in y_labels])
    
    plt.title('T = '+ str(real_year) + ' years')
    plt.grid(b=True, axis='x', which='major', color='black', linestyle='-',linewidth=0.2,alpha=0.5)
    plt.grid(b=True, axis='y', which='major', color='black', linestyle='-',linewidth=0.2,alpha=0.5)    
    
    figure = plt.gcf()
    figure.set_size_inches(4, 4)
    export_file = str(loop_idx) + '_' + str(real_year) + '_hist.png'
    plt.savefig(export_file,dpi=150,bbox_inches = "tight")
    
    plt.clf()
    plt.close('all')
    