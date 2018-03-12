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



posi = np.array([[1,2,3],[4,5,6]])
velo = 12*np.array([[1,2,3],[4,5,6]])


print posi
print velo

print np.cross(posi, velo)


x = np.array([[1,2,3], [4,5,6], [7, 8, 9]])
y = np.array([[7, 8, 9], [4,5,6], [1,2,3]])
print np.cross(x, y)
