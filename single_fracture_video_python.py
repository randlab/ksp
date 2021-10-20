#PYTHON 3 (THREE)
#C:\Users\A480325\Miniconda3\python.exe
#C:\Users\A480325\Miniconda3\python.exe -m pip install moviepy

import os
import moviepy.video.io.ImageSequenceClip
import re
from PIL import Image
from tqdm import tqdm

frac_dir     = 'C:/Daten/KARST_NAGRA/report/python_figure/figures_1K_single_fracture/1_fracs'
histo_dir    = 'C:/Daten/KARST_NAGRA/report/python_figure/figures_1K_single_fracture/2_histo'
flowEq_dir   = 'C:/Daten/KARST_NAGRA/report/python_figure/figures_1K_single_fracture/3_flowEq'
flowRt_dir   = 'C:/Daten/KARST_NAGRA/report/python_figure/figures_1K_single_fracture/4_flowRt'
merged_dir   = 'C:/Daten/KARST_NAGRA/report/python_figure/figures_1K_single_fracture/5_merge'

flowEq_files = [flowEq_dir+'/'+img for img in os.listdir(flowEq_dir) if img.endswith("png")]
flowEq_files.sort(key=lambda f: int(re.sub('\D', '', f)))

flowRt_files = [flowRt_dir+'/'+img for img in os.listdir(flowRt_dir) if img.endswith(".png")]
flowRt_files.sort(key=lambda f: int(re.sub('\D', '', f)))

frac_files = [frac_dir+'/'+img for img in os.listdir(frac_dir) if img.endswith(".png")]
frac_files.sort(key=lambda f: int(re.sub('\D', '', f)))

histo_files = [histo_dir+'/'+img for img in os.listdir(histo_dir) if img.endswith(".png")]
histo_files.sort(key=lambda f: int(re.sub('\D', '', f)))

#for pic in frac_files:
#    print (pic)

print ('Merging images')
for index, file in enumerate (tqdm (frac_files, unit = ' images') ):

    image1 = Image.open(file)                   #frac
    image2 = Image.open(histo_files[index])#    #histo
    image3 = Image.open(flowEq_files[index])    #flowEq
    image4 = Image.open(flowRt_files[index])    #flowRt

    image1 = image1.resize((1600,602))
    image3 = image3.resize((1600,602))
    
    image1_size = image1.size
    image2_size = image2.size
    image3_size = image3.size
    image4_size = image4.size
       
    #♣new_image = Image.new('RGB',(2187, 2*image1_size[1]), (250,250,250))
    new_image = Image.new('RGB',(2215, 2*image1_size[1]), (250,250,250))
    new_image.paste(image1,(0,0))
    new_image.paste(image2,(image1_size[0],0))
    new_image.paste(image3,(0,image1_size[1]))
    new_image.paste(image4,(image1_size[0]+23,image1_size[1]))
    save_name = os.path.join (merged_dir,str(index)+".png")
    new_image.save(save_name,"PNG")

merged_files = [merged_dir+'/'+img for img in os.listdir(merged_dir) if img.endswith(".png")]
merged_files.sort(key=lambda f: int(re.sub('\D', '', f)))
      
print ('Assembling clip')    
fps=1
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(merged_files, fps=fps)
clip.write_videofile('my_video.mp4')