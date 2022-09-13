import pandas as pd
testData = pd.read_csv('isingTestData0.txt', skiprows=1, header=None)
testData = pd.DataFrame(testData.values.reshape(-1, 9))
testData.to_csv('isingTestData.csv')

import scipy.linalg as la
import numpy as np
offdi = la.circulant([0, 1, 1])
I = np.eye(3)

A = np.kron(offdi,I) + np.kron(I,offdi)

pd.DataFrame(A.astype(int)).to_csv('isingTestData_trueAdjMat.csv')