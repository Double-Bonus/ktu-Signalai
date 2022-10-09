# student number = 10
# accord = Dm

# f1 = 0 # is not used

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.io.wavfile import write
from sklearn import preprocessing

#TODO: implement logger instead of if debug

class MusicLab:
    def __init__(self, debug = False):
        self.debug = debug

        f2 = 110
        f3 = 147
        f4 = 220
        f5 = 294
        f6 = 349
        self.notes = [f2, f3, f4, f5, f6]
        self.notesNames = ["A styga", "D styga", "G styga", "B styga", "e styga"]
        self.STRING_COUNT = 5
        self.saveDir = "out/"
        
        self.samplingRate = 44100
        self.t_s = 3
        self.N = np.empty(self.STRING_COUNT, dtype=int)

    def countDelays_N(self):
        for i, note in enumerate(self.notes):
            if self.debug:
                print(i)
                print(note)
            self.N[i] = round(self.samplingRate / note) # 1 Task
            
    def generateInput_X(self): # 2 Task
        x_final = []
        x_random = []
        
        for i in range(0, self.STRING_COUNT):
            x_random = np.random.uniform(0, 1, self.N[i])
            K_zeroCount = self.samplingRate*self.t_s - self.N[i]
            if self.debug:
                print(K_zeroCount)

            x_zeros = []
            x_zeros = np.zeros(K_zeroCount)

            x_final.append(np.concatenate([x_random, x_zeros]))
            if self.debug:
                print(x_final[i])
        return x_final

    def generateSound_Y(self, signal_x):
        y_final = []
        for i in range(0, self.STRING_COUNT):
            b = [1]
            a = np.concatenate([[1], np.zeros(self.N[i]-1), [-0.5, -0.5]])
            
            audioData = signal.lfilter(b, a, signal_x[i])
            audioScaled = preprocessing.minmax_scale(audioData, feature_range=(-1,1))

            y_final.append(audioScaled); # 3 Task

        if self.debug:
            print(np.shape(signal_x))
        # static assert ar len == (fd*ts=44100*3=132300)
        return y_final
    
    # 4. Listen to notes:
    # TODO: sampling rate no need as parameter????
    def saveNoteAsWav(self, noteData, samplingRate, filename):
        write(filename=filename, rate=samplingRate, data=noteData.astype(np.float32))
        
    def drawSignal(self, signal_y, title, show=False):
        tn = np.linspace(0, self.t_s, num=len(signal_y))
        plt.figure
        plt.plot(tn, signal_y, 'r-')
        plt.title('Signalas laiko srityje, ' + title)
        plt.xlabel('t, s')
        plt.ylabel('A')
        plt.grid(True)
        plt.savefig(self.saveDir + "amp_" + title)
        if show:
            plt.show()
        plt.close()

    # 5 Get FFT and draw spectrum
    def drawSpectrum(self, signal_y, title, show = False):
        nfft = len(signal_y)
        yf = np.fft.fft(signal_y)

        spectrum = np.abs(yf) / nfft
        spectrum_db = 20 * np.log10(spectrum/np.max(spectrum))
        if 0: # debug TODO: redo
            print(yf)
            print(type(yf))
            print(spectrum)
            print(spectrum_db)
            
        k = list(range(0, nfft))
        f_Hz = [i * (self.samplingRate/nfft) for i in k] 

        # fig = plt.figure()
        ax = plt.axes()
        ax.plot(f_Hz, spectrum_db)
        ax.set_xlim(0, self.samplingRate/2)
        ax.set_ylim(-80, 0)
        plt.title('Signalas dažnių srityje, ' + title)
        plt.xlabel('f, Hz')
        plt.ylabel('S, db')
        plt.grid(True)
        plt.savefig(self.saveDir + "spec_" + title)
        if show:
            plt.show()
        plt.close()
        
    # 3.1.2 Simulate accord
    def generateAccord(self, allNotes):
        delay_ms = 75
        second_ms = 1000
        n_delay =  (delay_ms * self.samplingRate / second_ms)
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
    
    # 3.2.1
    def nonLinearDistortion(self, signal):
        def satlins(n):
            if (n <= -1):
                return -1
            if (-1 <= n <= 1):
                return n
            if (1 <= n):
                return 1
        sig_after = [satlins(i) for i in signal]
        return sig_after
    
    """ 3.2.2 Modeling the reverberation effect
        # filter coefficients
        # b = [1]
        # a = [1; 0...0(N); -(K)]
    """
    def addReverb(self, signalIn, N_ms, K_coef):
        second_ms = 1000
        n_delay =  (N_ms * self.samplingRate / second_ms)
        b = [1]
        a = np.concatenate([[1], np.zeros(int(n_delay-1)), [-K_coef]])
        reverbedSignal = signal.lfilter(b, a, signalIn)
        return reverbedSignal

    
""" 4. a) 
Overdrive is an effect where the amplitude of the input signal undergoes a non-linear amplification. The threshold determines how much of 
the signal undergoes the nonlinear amplification curve (lower threshold implies more of the signal is captured in the curve).
"""
def overdrive(signal):
    def overLambda(x):
        # Remove a result of floating-point approximation for python floats
        if x > 1 and x < 1.000000000000004:
            x = 1
        elif x < -1 and x > -1.000000000000004:
            x = -1
        
        if (0 <= np.abs(x) and np.abs(x) < (1/3)):
            return x * 2 * np.abs(x)
        elif ((1/3) <= np.abs(x) and np.abs(x) < (2/3)):
            return x * (3-(2-(3*np.abs(x))**2 )) / 3
        elif ((2/3) <= np.abs(x) and np.abs(x) <= 1):
            return x
        else:
            print("Incorrect value!!!")
            print(x)

    sig_after = [overLambda(i) for i in signal if overLambda(i) is not None]
    return sig_after

def applyFuzz(signal, a = 50):
    fuzz = lambda x: x * (1 - np.exp(-a * np.abs(x)))
    sig_after = fuzz(signal)
    return sig_after
            
############### main

plt.rcParams.update({'font.family': "Times New Roman"})
plt.rcParams.update({'font.size': 10})

musicObj = MusicLab(True)
musicObj.countDelays_N()
print(f"Signal delays {musicObj.N}")

signals_X = musicObj.generateInput_X()
sounds_Y = musicObj.generateSound_Y(signals_X)

for i, audioData in enumerate(sounds_Y):
    musicObj.saveNoteAsWav(audioData, musicObj.samplingRate, f"Note{i}.wav")

if 1:
    for i, y_sig in enumerate(sounds_Y):
        musicObj.drawSignal(y_sig, musicObj.notesNames[i])
        musicObj.drawSpectrum(y_sig, musicObj.notesNames[i])

if 1:
    accordSignal = musicObj.generateAccord(sounds_Y)
    musicObj.drawSignal(accordSignal, "Dm akordas", True)
    musicObj.drawSpectrum(accordSignal, "Dm akordas", True)
    musicObj.saveNoteAsWav(accordSignal, musicObj.samplingRate, "accord.wav")


# Make wrapper for these 5 functions??
if 0:
    K = 30
    distortedAccord = musicObj.nonLinearDistortion(accordSignal)
    distortedAccord = np.multiply(distortedAccord, K)

    musicObj.drawSignal(distortedAccord)
    musicObj.drawSpectrum(distortedAccord)
    musicObj.saveNoteAsWav(distortedAccord, musicObj.samplingRate, f"DistAccord{K}.wav")

# Analyze how the sound of the chord and its temporal and frequency characteristics change when K = 5 and K = 50.
if 0:
    K = 5
    distortedAccord_5 = musicObj.nonLinearDistortion(accordSignal)
    distortedAccord_5 = np.multiply(distortedAccord_5, K)
    musicObj.drawSignal(distortedAccord_5)
    musicObj.drawSpectrum(distortedAccord_5)
    musicObj.saveNoteAsWav(distortedAccord_5, musicObj.samplingRate, f"DistAccord{K}.wav")

    K = 50
    distortedAccord_50 = musicObj.nonLinearDistortion(accordSignal)
    distortedAccord_50 = np.multiply(distortedAccord_50, K)
    musicObj.drawSignal(distortedAccord_50)
    musicObj.drawSpectrum(distortedAccord_50)
    musicObj.saveNoteAsWav(distortedAccord_50, musicObj.samplingRate, f"DistAccord{K}.wav")


# 3.2.2 task
if 0:
    N_ms = 200
    K_reverb = 0.5

    accordSignal = musicObj.generateAccord(sounds_Y)
    accordSignal_rev = musicObj.addReverb(accordSignal, N_ms, K_reverb)

    musicObj.drawSignal(accordSignal_rev, musicObj.samplingRate)
    musicObj.drawSpectrum(accordSignal_rev)
    musicObj.saveNoteAsWav(accordSignal_rev, musicObj.samplingRate, f"Reverb{K_reverb}.wav")


# 4 extra task
if 0:
    newSounds_Y = musicObj.generateSound_Y(signals_X)
    overSounds_Y = []
    for sound in newSounds_Y:
        print("tests")
        overSounds_Y.append(overdrive(sound))
        
    accordSignal_o = musicObj.generateAccord(sounds_Y)
    musicObj.drawSignal(accordSignal_o)
    musicObj.drawSpectrum(accordSignal_o)
    
    newAccordSignal_o = musicObj.generateAccord(overSounds_Y)
    musicObj.drawSignal(newAccordSignal_o) # TODO: check if thats really zero
    musicObj.drawSpectrum(newAccordSignal_o)
    musicObj.saveNoteAsWav(newAccordSignal_o, musicObj.samplingRate, f"OverNew.wav")
        

    # accordSignal_o = musicObj.generateAccord(sounds_Y)
    # musicObj.drawSignal(accordSignal_o)
    # musicObj.drawSpectrum(accordSignal_o)
    # accordSignal_over = overdrive(accordSignal_o)

    # musicObj.drawSignal(accordSignal_over) # TODO: check if thats really zero
    # musicObj.drawSpectrum(accordSignal_over)
    # musicObj.saveNoteAsWav(accordSignal_over, musicObj.samplingRate, f"Over.wav")


if 0:
    accordSignal_f = musicObj.generateAccord(sounds_Y)
    musicObj.drawSignal(accordSignal_f)
    musicObj.drawSpectrum(accordSignal_f)
    
    accordSignal_fuzz = applyFuzz(accordSignal_f, a=15) # or to notes?!?!?!?

    musicObj.drawSignal(accordSignal_fuzz)
    musicObj.drawSpectrum(accordSignal_fuzz)
    musicObj.saveNoteAsWav(accordSignal_fuzz, musicObj.samplingRate, f"Fuzz.wav")
