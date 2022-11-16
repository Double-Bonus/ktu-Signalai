import scipy.io
import matplotlib.pyplot as plt
import numpy as np


matData = scipy.io.loadmat('3-lab/lab3_signalai.mat')

print(matData)
print(type(matData))
variklioSig = matData.get('variklioSig')
kabinosSig = matData.get('kabinosSig')
pilotoSig = matData.get('pilotoSig')

print(variklioSig)
print(type(variklioSig))
print(variklioSig.shape)

print(kabinosSig)
print(type(kabinosSig))

# plt.plot(kabinosSig, color="red")
# plt.show()


a = np.zeros(176000)
print(type(a))

plt.plot(variklioSig)
plt.show()

