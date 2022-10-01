# student number = 10
# accord = Dm

# f1 = 0 # is not used

import numpy as np
from scipy import signal
from scipy.io.wavfile import write
import matplotlib.pyplot as plt
from sklearn import preprocessing

debug = True

f2 = 110
f3 = 147
f4 = 220
f5 = 294
f6 = 349

Fd  = 44100
t_s = 3

notes = [f2, f3, f4, f5, f6]


N=[]
X_final = []
y_final = []
for i, note in enumerate(notes):
    if debug:
        print(i)
        print(note)
    N.append(round(Fd / note)) # 1 Task

    if debug:
        print(N)

    x_random = []
    x_random = np.random.uniform(0, 1, N[i])
    if debug:
        print(x_random)

    K_zeroCount = Fd*t_s - N[i]
    if debug:
        print(K_zeroCount)

    x_zeros = []
    x_zeros = np.zeros(K_zeroCount)


    X_final.append(np.concatenate([x_random, x_zeros])) # 2 Task
    if debug:
        print(X_final[i])

    b = [1]
    a = np.concatenate([[1], np.zeros(N[i]-1), [-0.5, -0.5]])
    
    y_final.append(signal.lfilter(b, a, X_final[i])); # 3 Task

print(np.shape(X_final))
# static assert ar len == (fd*ts=44100*3=132300)


tn = np.linspace(0, t_s, num=(Fd*t_s))
if 0:
    plt.figure
    plt.plot(tn, X_final[0], 'b', alpha=0.75)
    plt.plot(tn, y_final[0], 'r--')
    # plt.plot(tn, y_final[0], 'r--', t, z2, 'r', t, y, 'k')
    # plt.legend(('noisy signal', 'lfilter, once', 'lfilter, twice',
    #             'filtfilt'), loc='best')
    plt.grid(True)
    plt.show()


# 4. Listen to notes:
def saveNoteAsWav(noteData, samplingRate, filename):
    audioScaled = preprocessing.minmax_scale(noteData, feature_range=(-1,1))
    write(filename=filename, rate=Fd, data=audioScaled.astype(np.float32))

for i, audioData in enumerate(y_final):
    saveNoteAsWav(audioData, Fd, f"Note{i}.wav")
    

# 5 Get FFT and see spectrum


noteNr = 2

nfft = len(y_final[noteNr])
yf = np.fft.fft(y_final[noteNr])

spectrum = np.abs(yf) / nfft
spectrum_db = 20 * np.log10(spectrum/np.max(spectrum))

if 0: # debug
    print(yf)
    print(type(yf))
    print(spectrum)
    print(spectrum_db)
    
k = list(range(0, nfft))

print("5 Debug")
print(Fd)
print(nfft)

f_Hz = [i * (Fd/nfft) for i in k]



if 1:    
    # fig, ax= plt.figure()
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(f_Hz, spectrum_db)
    ax.set_xlim(0, Fd/2)
    ax.set_ylim(-80, 0)
    # plt.plot(tn, y_final[0], 'r--')
    # plt.plot(tn, y_final[0], 'r--', t, z2, 'r', t, y, 'k')
    # plt.legend(('noisy signal', 'lfilter, once', 'lfilter, twice',
    #             'filtfilt'), loc='best')
    ax.grid(True)
    plt.show()






# nfft = len(X_final[0])
# S = abs()