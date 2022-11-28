import scipy.io
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile



class Adapt:
    def __init__(self, debug = False):
        self.debug = debug
        self.time_s = 22
        self.Fs_Hz = 8000
        self.saveDir = "out/"
         

    def drawSignal(self, signal_y, title, show=False):
            tn = np.linspace(0, self.time_s, num=len(signal_y))
            plt.figure
            plt.plot(tn, signal_y, 'r-')
            plt.title('Signalas laiko srityje, ' + title)
            plt.xlabel('t, s')
            plt.ylabel('A')
            plt.grid(True)
            plt.savefig(self.saveDir + "amp_" + title,  bbox_inches='tight',
                        pad_inches=0)
            if show:
                plt.show()
            plt.close()

    # 5 FFT and draw spectrum
    def drawSpectrum(self, signal_y, title, show= False):
        nfft = len(signal_y)
        yf = np.fft.fft(signal_y)

        spectrum = np.abs(yf) / nfft
        spectrum_db = 20 * np.log10(spectrum/np.max(spectrum))

        k = list(range(0, nfft))
        f_Hz = [i * (self.Fs_Hz/nfft) for i in k]

        ax = plt.axes()
        ax.plot(f_Hz, spectrum_db)
        ax.set_xlim(0, self.Fs_Hz/2)
        ax.set_ylim(-80, 0)
        plt.title('Signalas dažnių srityje, ' + title)
        plt.xlabel('f, Hz')
        plt.ylabel('S, db')
        plt.grid(True)
        plt.savefig(self.saveDir + "spec_" + title,  bbox_inches = 'tight',
                    pad_inches= 0)
        if show:
            plt.show()
        plt.close()
    
    def saveSignalAsWav(self, noteData, filename):
        wavfile.write(filename=filename, rate=self.Fs_Hz, data=noteData.astype(np.float32))
   
 
# 3.4.1 MVK algoritmas
def calculate_MVK(inNoise, inCombine, filterOrder=20, step = 0.1):

    if len(inCombine) != len(inNoise):
        print("Signals are not the same length")
        exit -1
        
    w   = np.transpose(np.zeros((1, filterOrder)))
    x_a = np.transpose(np.zeros((1, filterOrder)))
    
    s_iv = np.zeros((len(inNoise),))

    for n in range(len(inNoise)):
        x_a = np.roll(x_a, 1)
        x_a[0] = inNoise[n]
        x_iv = np.matmul(np.transpose(w), x_a)
        s_iv[n] = inCombine[n] - x_iv[0]
        w = w + 2*step*s_iv[n]*x_a
    return s_iv

def calculate_MSE(real, prediction):
    mse = (np.square(real - prediction)).mean()
    return mse

matData = scipy.io.loadmat('signalai/lab3_signalai.mat')

print(matData)
print(type(matData))
variklioSig = matData.get('variklioSig')[0]
kabinosSig = matData.get('kabinosSig')[0]
pilotoSig = matData.get('pilotoSig')[0]

print(type(kabinosSig))

plt.rcParams.update({'font.family': "Times New Roman"})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.figsize': (16, 6)}) # horizontally longer figure 

filterObj = Adapt()

if 0:
    # 3.3.1
    filterObj.drawSignal(variklioSig, "variklis")
    filterObj.drawSignal(kabinosSig, "kabinos")
    filterObj.drawSignal(pilotoSig, "piloto")

    filterObj.drawSpectrum(variklioSig, "Variklis")
    filterObj.drawSpectrum(kabinosSig, "kabinos")
    filterObj.drawSpectrum(pilotoSig, "piloto")


    # 3.3.2
    filterObj.saveSignalAsWav(variklioSig, f"out/variklis.wav")
    filterObj.saveSignalAsWav(kabinosSig, f"out/kabina.wav")
    filterObj.saveSignalAsWav(pilotoSig, f"out/pilotas.wav")


if 0:
    # 3.4 
    sig_afterMVK = calculate_MVK(variklioSig, kabinosSig)

    filterObj.drawSignal(sig_afterMVK, "sig_afterMVK")
    filterObj.drawSpectrum(sig_afterMVK, "sig_afterMVK")
    filterObj.saveSignalAsWav(sig_afterMVK, f"out/sigMVK.wav")



if 0:
    # 3.5.6
    M_filterOrder = 20
    mu_step = 0.1
    M_array = [M_filterOrder, M_filterOrder*2, M_filterOrder//2]
    mu_string = str(mu_step)
    mu_string = mu_string.replace('.', ',')  # Convert for saving
    print("Analysis with different filter order")
    for M in M_array:
        sig_MVK = calculate_MVK(variklioSig, kabinosSig, M, mu_step)
        filterObj.drawSignal(sig_MVK, f"piloto signalas, kai M = {M}, mu = {mu_string}")
        filterObj.drawSpectrum(sig_MVK, f"piloto signalas, kai M = {M}, mu = {mu_string}")
        filterObj.saveSignalAsWav(sig_MVK, f"piloto signalas, kai M = {M}, mu = {mu_string}")
    


# 3.5.7
# Filtro koeficientų skaičių 𝑀 galite keisti logaritminiu masteliu: 10,
# 20, 50, 100. Adaptacijos žingsnio 𝜇 vertę galite keisti dėsniu 0.001, 0.005, 0.01, 0.05, 0.1.

sig_MVK = calculate_MVK(variklioSig, kabinosSig)

print(calculate_MSE(pilotoSig, sig_MVK))

M = [10, 15, 20, 25, 50, 100] # there are better ways
# M = [10, 15, 20] # there are better ways

u_mu = [0.001, 0.005, 0.01, 0.05, 0.1]
# u_mu = [0.001, 0.01, 0.1]

MSE_array = np.zeros((len(M), len(u_mu)))

for i in range(len(M)):
    for j in range(len(u_mu)):
        MSE_array[i, j] = calculate_MSE(pilotoSig, calculate_MVK(variklioSig, kabinosSig, M[i], u_mu[j]))

# find min - i and j

         
plt.figure
plt.plot(M, MSE_array) # check markers
plt.title('MSE nuo M ir mu')
plt.xlabel('M')
plt.ylabel('MSE')
plt.grid(True)
# plt.savefig(self.saveDir + "amp_" + title,  bbox_inches='tight', pad_inches=0)
# add legend
plt.show()


# plot MSE_array over time



