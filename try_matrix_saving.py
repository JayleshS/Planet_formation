import numpy as np

matrix = np.array([1])
save_matrix = np.array([])
save_matrix_zeros = np.zeros(2)
save_list = []
i= 1

while i < 5.:

  matrix += 1

  save_list.append(matrix)
  i += 1

print 'save_list: ',save_list

print type(matrix)
