import matplotlib.pyplot as plt
import numpy as np


# parameters
w = 17
halfw = w // 2           # floor(w/2)
h = 25
halfh = h // 2           # floor(h/2)

# generate the test image
inpaintI = np.ones((100, 100), dtype=float)

cols_start_mat = 50 - halfw + 1   
cols_end_mat   = 50 - halfw + w  
rows_slice1 = slice(11 - 1, 90)                 
cols_slice1 = slice(cols_start_mat - 1, cols_end_mat)  

inpaintI[rows_slice1, cols_slice1] = 0.0

# generate mask
mask = np.zeros((100, 100), dtype=bool)

rows_start_mat = 50 - halfh + 1   
rows_end_mat   = 50 - halfh + h  
rows_slice2 = slice(rows_start_mat - 1, rows_end_mat)  
cols_slice2 = slice(21 - 1, 80)                        

mask[rows_slice2, cols_slice2] = True

# Apply mask
inpaintI[mask] = 0.5

# Display
plt.figure()
plt.imshow(inpaintI, cmap='gray', vmin=0, vmax=1)
plt.axis('off')
plt.title('inpaintI')

plt.figure()
plt.imshow(mask, cmap='gray')
plt.axis('off')
plt.title('mask')

plt.show()
