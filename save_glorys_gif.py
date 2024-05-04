#!/usr/bin/env python
# coding: utf-8

# In[3]:


import glob
from PIL import Image
from natsort import natsorted


# In[ ]:





# In[5]:

depth = 1500
slow = 500
fast = 200

fnames = glob.glob(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/*_{depth}m.png")   
files = natsorted(fnames)
print(files)


# In[ ]:


# Create the frames
frames = []
imgs = files[:]
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)
 
# Save into a GIF file that loops forever
frames[0].save(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/1993_2019_{depth}m_fast.gif", format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=fast, loop=0)


# In[ ]:




