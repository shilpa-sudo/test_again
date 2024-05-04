#!/usr/bin/env python
# coding: utf-8

# In[3]:


import glob
from PIL import Image
from natsort import natsorted


# In[ ]:


year = 1993
month = 1
slow = 500
fast = 200


# In[5]:


#fnames = glob.glob(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/*_0m.png")  for dates animation
fnames = glob.glob(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/1993-01_*")  #for depths animation
files = natsorted(fnames)
print(files)


# In[7]:


# Create the frames
frames = []
imgs = files[:]
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)
 
# Save into a GIF file that loops forever
#frames[0].save(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/1993_2019_0m.gif", format='GIF',
frames[0].save(f"/home/shilpa/glory_mat_analysis/vertical_anal_gl_quarter/{year}_{month}_alldepths.gif", format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=slow, loop=0)


# In[ ]:




