#%% get the data from vrlab  quiver animation
import numpy as np
import glob as glob

a1 = np.arange(1,61,1)
b1 = np.arange(1,61,1)

a2 = np.arange(1,61,1).reshape(12,5)
b2 = np.arange(1,61,1).reshape(12,5)

a3 = np.arange(1,61,1).reshape(3,4,5)
b3 = np.arange(1,61,1).reshape(3,4,5)
# %%
print(a1.shape, b1.shape)
print(np.stack((a1, b1)).shape)
print(np.vstack((a1, b1)).shape)
print(np.hstack((a1, b1)).shape)
print(np.dstack((a1, b1)).shape)

#%%
print(a2.shape, b2.shape)
print(np.stack((a2, b2)).shape)
print(np.vstack((a2, b2)).shape)
print(np.hstack((a2, b2)).shape)
print(np.dstack((a2, b2)).shape)
#%%
print(a3.shape, b3.shape)
print(np.stack((a3, b3)).shape)
print(np.vstack((a3, b3)).shape)
print(np.hstack((a3, b3)).shape)
print(np.dstack((a3, b3)).shape)

#%% we can use reshape to change the dimension
np.shape(np.reshape(a2, (-1,3,4)))

# %%
