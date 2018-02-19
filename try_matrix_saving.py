import numpy as np

matrix = np.array([1])
save_matrix = np.array([])
save_matrix_zeros = np.zeros(2)
save_list = []
i= 1

while i < 5.:

  matrix += 1
  np.append(save_matrix, matrix)
  np.append(save_matrix_zeros, matrix)
  save_list.append(list(matrix)))
  i += 1

print 'matrix ',matrix
print 'save_matrix: ', save_matrix
print 'save_matrix_zeros: ', save_matrix_zeros
print 'save_list: ',save_list
