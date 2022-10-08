# student number = 10
# accord = Dm

# f1 = 0 # is not used

import sys
import numpy as np
from scipy import signal
from scipy.io.wavfile import write
import matplotlib.pyplot as plt
from sklearn import preprocessing

debug = False

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
    
    
    audioData = signal.lfilter(b, a, X_final[i])
    audioScaled = preprocessing.minmax_scale(audioData, feature_range=(-1,1))
    
    y_final.append(audioScaled); # 3 Task

if debug:
    print(np.shape(X_final))
# static assert ar len == (fd*ts=44100*3=132300)

# 4. Listen to notes:
def saveNoteAsWav(noteData, samplingRate, filename):
    write(filename=filename, rate=samplingRate, data=noteData.astype(np.float32))

for i, audioData in enumerate(y_final):
    saveNoteAsWav(audioData, Fd, f"Note{i}.wav")
    

# 5 Get FFT and draw spectrum
def drawSignal(signal_y, t_s, Fd):
    tn = np.linspace(0, t_s, num=(Fd*t_s))
    plt.figure
    # ax = plt.axes()
    # ax.margins(0.2, 0.2)
    plt.plot(tn, signal_y, 'r--') #  alpha=0.75)
    # plt.plot(tn, y_final[0], 'r--', t, z2, 'r', t, y, 'k')
    # plt.legend(('noisy signal', 'lfilter, once', 'lfilter, twice',
    #             'filtfilt'), loc='best')
    plt.grid(True)
    plt.show()

def drawSpectrum(signal_y):
    nfft = len(signal_y)
    yf = np.fft.fft(signal_y)

    spectrum = np.abs(yf) / nfft
    spectrum_db = 20 * np.log10(spectrum/np.max(spectrum))
    if 0: # debug
        print(yf)
        print(type(yf))
        print(spectrum)
        print(spectrum_db)
        
    k = list(range(0, nfft))

    if 0: # debug
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
        # plt.legend(('noisy signal', 'lfilter, once', 'lfilter, twice',
        #             'filtfilt'), loc='best')
        ax.grid(True)
    plt.show()


# for i, y_sig in enumerate(y_final):
#     drawSignal(y_sig, t_s, Fd)
#     drawSpectrum(y_sig)



# 3.1.2 Simulate accord
def generateAccord(allNotes, samplingRate):
    delay_ms = 75
    second_ms = 1000
    n_delay =  (delay_ms * samplingRate / second_ms)
    # delay each note for some time (delay time: first t, second t*2, t*3, .... t*N)
    
    notesWithDelay=[]
    for i, note in enumerate(allNotes):
        if i == 0: # do not delay first note
            notesWithDelay.append(np.concatenate([note]))
        else:
            count = int(i * n_delay)
            temp_a = note[:-count]
            delayZeros = np.zeros(count).astype(np.float64)
            notesWithDelay.append(np.concatenate([delayZeros, temp_a]))
    
    accord = np.zeros(len(allNotes[0])).astype(np.float64)
    for note in notesWithDelay:
        accord = accord + note
    
    return accord


accordSignal = generateAccord(y_final, Fd)
# drawSignal(accordSignal, t_s, Fd)
# drawSpectrum(accordSignal)
saveNoteAsWav(accordSignal, Fd, "accord.wav")


# 3.2.1
def nonLinearDistortion(signal):
    def satlins(n):
        if (n <= -1):
            return -1
        if (-1 <= n <= 1):
            return n
        if (1 <= n):
            return 1
    sig_after = [satlins(i) for i in signal]
    return sig_after
    
K = 30
distortedAccord = nonLinearDistortion(accordSignal)
distortedAccord = np.multiply(distortedAccord, K) # is it good?

# drawSignal(distortedAccord, t_s, Fd)
# drawSpectrum(distortedAccord)
# saveNoteAsWav(distortedAccord, Fd, f"DistAccord{K}.wav")


# Aptarkite kaip keičiasi akordo skambesys ir jo laikinės bei dažninės charakteristikos, kai K = 5 ir K = 50.
if 0:
    K = 5
    distortedAccord_5 = nonLinearDistortion(accordSignal)
    distortedAccord_5 = np.multiply(distortedAccord_5, K) # is it good?

    drawSignal(distortedAccord_5, t_s, Fd)
    drawSpectrum(distortedAccord_5)
    saveNoteAsWav(distortedAccord_5, Fd, f"DistAccord{K}.wav")

    K = 50
    distortedAccord_50 = nonLinearDistortion(accordSignal)
    distortedAccord_50 = np.multiply(distortedAccord_50, K) # is it good?

    drawSignal(distortedAccord_50, t_s, Fd)
    drawSpectrum(distortedAccord_50)
    saveNoteAsWav(distortedAccord_50, Fd, f"DistAccord{K}.wav")


# 3.2.2 Reverberacijos efekto modeliavimas
N_ms = 200
K_reverb = 1

# filtro koeficientai:
# b = [1]
# a = [1; 0...0(N);  0.4(K)]

def addReverb(signalIn, samplingRate, N_ms, K_coef):
    second_ms = 1000
    n_delay =  (N_ms * samplingRate / second_ms)
    b = [1]
    a = np.concatenate([[1], np.zeros(int(n_delay-1)), [-K_coef]])
    reverbedSignal = signal.lfilter(b, a, signalIn)
    return reverbedSignal

# print("reverbedSignal")
# accordSignal = generateAccord(y_final, Fd)
# accordSignal_rev = addReverb(accordSignal, Fd, N_ms, K_reverb)

# drawSignal(accordSignal_rev, t_s, Fd)
# drawSpectrum(accordSignal_rev)
# saveNoteAsWav(accordSignal_rev, Fd, f"Reverb{K_reverb}.wav")


# 4 extra task
if 0:
    """ 
    Overdrive is an effect where the amplitude of the input signal undergoes a non-linear amplification. The threshold determines how much of 
    the signal undergoes the nonlinear amplification curve (lower threshold implies more of the signal is captured in the curve).
    """
    def overdrive(signal):
        def overLambda(x):
            if (0 <= np.abs(x) and np.abs(x) < (1/3)):
                return x * 2 * np.abs(x)
            elif ((1/3) <= np.abs(x) and np.abs(x) < (2/3)):
                return x * (3-(2-(3*np.abs(x))**2 )) / 3
            elif ((2/3) <= np.abs(x) and np.abs(x) < 1):
                return x
        sig_after = [overLambda(i) for i in signal if overLambda(i) is not None]
        return sig_after

    accordSignal_o = generateAccord(y_final, Fd)
    drawSignal(accordSignal_o, t_s, Fd)
    drawSpectrum(accordSignal_o)
    accordSignal_over = overdrive(accordSignal_o)

    drawSignal(accordSignal_over, t_s, Fd)
    drawSpectrum(accordSignal_over)
    saveNoteAsWav(accordSignal_over, Fd, f"Over.wav")


if 1:
    def applyFuzz(signal, a = 50):
        def fuzzLambda(x):
            return x * (1 - np.exp( -a* np.abs(x)))
        fuzz = lambda x: x * (1 - np.exp( -a* np.abs(x)))
        # sig_after = [fuzzLambda(i) for i in signal]
        sig_after = fuzz(signal)
        return sig_after

    accordSignal_f = generateAccord(y_final, Fd)
    accordSignal_fuzz = applyFuzz(accordSignal_f)

    drawSignal(accordSignal_fuzz, t_s, Fd)
    drawSpectrum(accordSignal_fuzz)
    saveNoteAsWav(accordSignal_fuzz, Fd, f"Fuzz.wav")
