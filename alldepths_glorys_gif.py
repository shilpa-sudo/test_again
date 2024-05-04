#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import glob
from PIL import Image
from natsort import natsorted
import numpy as np


# In[ ]:


slow = 500
fast = 200


# In[ ]:


for k in range(1993,2020):
    for i in range(1,10):
        y = f'{k}-0{i}'
        
        fnames = glob.glob(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/{y}_*")  #for depths animation
        files = natsorted(fnames)
        
        
        # Create the frames
        frames = []
        imgs = files[:]
        for i in imgs:
            new_frame = Image.open(i)
            frames.append(new_frame)

        # Save into a GIF file that loops forever
        #frames[0].save(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/1993_2019_0m.gif", format='GIF',
        frames[0].save(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/{y}_alldepths.gif", format='GIF',
                       append_images=frames[1:],
                       save_all=True,
                       duration=slow, loop=0)


# In[ ]:


for k in range(1993,2020):
    for i in range(10,13):
        x = f'{k}-{i}'
        fnames = glob.glob(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/{x}_*")  #for depths animation
        files = natsorted(fnames)
        
        
        # Create the frames
        frames = []
        imgs = files[:]
        for i in imgs:
            new_frame = Image.open(i)
            frames.append(new_frame)

        # Save into a GIF file that loops forever
        #frames[0].save(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/1993_2019_0m.gif", format='GIF',
        frames[0].save(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/{x}_alldepths.gif", format='GIF',
                       append_images=frames[1:],
                       save_all=True,
                       duration=slow, loop=0)


# In[ ]:





# In[ ]:





# In[ ]:




