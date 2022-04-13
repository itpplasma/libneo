"""
  matrixutil.py
  assortment of useful routines on 2D matrices:
  - Bresenham line algorithm
  - floodfill an area defined by a boundary
  - 2D Thomas algorithm
"""

import numpy as np


# ..................................................................

def bresenham_line(x0, y0, x1, y1):
  """
    Bresenham line algorithm.
    Draw a line from [x0,y0] to [x1,y1].
    x0, y0, x1, y1 are integers, and indices to a 2D array (which is irrelevant here)
    Result is a list of tuples with [x,y] coordinates (indices) of the line produced.
  """
  steep = abs(y1 - y0) > abs(x1 - x0)
  if steep:
    x0, y0 = y0, x0  
    x1, y1 = y1, x1
  swapped = False
  if x0 > x1:
    swapped = True
    x0, x1 = x1, x0
    y0, y1 = y1, y0
  if y0 < y1: 
    ystep = 1
  else:
    ystep = -1
  deltax = x1 - x0
  deltay = abs(y1 - y0)
  error = -deltax / 2.0
  y = y0
  line = []    
  for x in range(x0, x1 + 1):
    if steep:
        line.append((y,x))
    else:
        line.append((x,y))
    error = error + deltay
    if error > 0:
        y = y + ystep
        error = error - deltax
  if swapped:
    line.reverse
  return line

# ..................................................................

def floodfill_recursive(matrix, x, y, v):        # caution highly recursive!
  """
    Flood fill an area in 2D 'matrix' with value v (v must not be 0).
    Start with indices [x,y].
    Touch only regions where matrix=0, i.e. the boundary is defined
    where matrix<>0.
    Caution: this routine is highly recursive, i.e. needs stack size
  """
  if matrix[x][y]==0:
    matrix[x][y] = v
    if x > 0:
      floodfill(matrix, x-1, y, v)
    if x < matrix.shape[0] - 1:
      floodfill(matrix, x+1, y, v)
    if y > 0:
      floodfill(matrix, x, y-1, v)
    if y < matrix.shape[1] - 1:
      floodfill(matrix, x, y+1, v)

      

def floodfill (data, start_coords, fill_value):
    """
    Flood fill algorithm
    
    Parameters
    ----------
    data : (M, N) ndarray of uint8 type
        Image with flood to be filled. Modified inplace.
    start_coords : tuple
        Length-2 tuple of ints defining (row, col) start coordinates.
    fill_value : int
        Value the flooded area will take after the fill.
        
    Returns
    -------
    None, ``data`` is modified inplace.
    """
    xsize, ysize = data.shape
    orig_value = data[start_coords[0], start_coords[1]]
    
    stack = set(((start_coords[0], start_coords[1]),))
    if fill_value == orig_value:
      raise ValueError("Filling region with same value "
                       "already present is unsupported. "
                       "Did you already fill this region?")
    while stack:
        x, y = stack.pop()

        if data[x, y] == orig_value:
            data[x, y] = fill_value
            if x > 0:
                stack.add((x - 1, y))
            if x < (xsize - 1):
                stack.add((x + 1, y))
            if y > 0:
                stack.add((x, y - 1))
            if y < (ysize - 1):
                stack.add((x, y + 1))

      
# ..................................................................
      
def tridag_multi (a, b, c, d):
  """
    Solve a set of tridiagonal systems.
  """
  nf = a.shape[0]     # number of equations
  ac, bc, cc, dc = map(np.array, (a, b, c, d))
  for it in xrange(1, nf):
    mc = ac[it,:] / bc[it-1,:]
    bc[it,:] = bc[it,:] - mc*cc[it-1,:] 
    dc[it,:] = dc[it,:] - mc*dc[it-1,:]
  xc = ac
  xc[-1,:] = dc[-1,:] / bc[-1,:]
  for il in xrange(nf-2, -1, -1):
    xc[il,:] = (dc[il,:]-cc[il,:]*xc[il+1,:])/bc[il,:]
  del bc, cc, dc  # delete variables from memory
  return xc


