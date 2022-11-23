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
    


matData = scipy.io.loadmat('signalai/lab3_signalai.mat')

print(matData)
print(type(matData))
variklioSig = matData.get('variklioSig')[0]
kabinosSig = matData.get('kabinosSig')[0]
pilotoSig = matData.get('pilotoSig')[0]


plt.rcParams.update({'font.family': "Times New Roman"})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.figsize': (16, 6)}) # horizontally longer figure 

filterObj = Adapt()

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




# 3.4.1 MVK algoritmas
filterOrder = 20
step = 0.1
# w = np.zeros(filterOrder, dtype=) 
w  = np.zeros(filterOrder) # Transpose
xa = np.zeros(filterOrder) 




# circshift
# >>> np.roll(a,2)
# array([8, 9, 0, 1, 2, 3, 4, 5, 6, 7])
# >>> np.roll(a,-2)
# array([2, 3, 4, 5, 6, 7, 8, 9, 0, 1])
