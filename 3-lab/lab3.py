import scipy.io
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from sklearn import preprocessing

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
            plt.ylabel('Amplitudė')
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

#3.6
def calculate_NMVK(inNoise, inCombine, filterOrder=20, step = 0.1):
    
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
        
        w = w + (step / (np.matmul(np.transpose(x_a), x_a))) * s_iv[n] * x_a # only this differs from MVK
    return s_iv

def calculate_RMK(inNoise, inCombine, filterOrder=15, Lam=0.99):
    
    if len(inCombine) != len(inNoise):
        print("Signals are not the same length")
        exit -1
    
    delta = 0.01
    I = np.eye(filterOrder)
    P = (delta**-1) * I
    
    w   = np.transpose(np.zeros((1, filterOrder)))
    x_a = np.transpose(np.zeros((1, filterOrder)))
    
    s_iv = np.zeros((len(inNoise),))
    
    for n in range(len(inNoise)): 
        x_a = np.roll(x_a, 1)
        x_a[0] = inNoise[n]
        
        v = np.matmul(P, x_a)
        u = np.matmul(np.transpose(P), v)
        
        v_norm = np.linalg.norm(v)
        k = 1 / (Lam + v_norm**2 + np.sqrt(Lam)*np.sqrt(Lam+v_norm**2))
        
        P = (P - k * np.matmul(v, np.transpose(u))) / np.sqrt(Lam)   
        
        x_iv = np.matmul(np.transpose(w), x_a)
        s_iv[n] = inCombine[n] - x_iv[0]
        if np.isnan(s_iv[n]):
            print(s_iv[n])
  
        w = w + ((s_iv[n] * u) / (Lam + v_norm**2))
    
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

# # Normalize in 1/-1 range
# variklioSig = preprocessing.minmax_scale(variklioSig, feature_range=(-1,1))
# kabinosSig = preprocessing.minmax_scale(kabinosSig, feature_range=(-1,1))
# pilotoSig = preprocessing.minmax_scale(pilotoSig, feature_range=(-1,1))

print("Loaded signals")
print(type(kabinosSig))

plt.rcParams.update({'font.family': "Times New Roman"})
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.figsize': (16, 6)}) # horizontally longer figure 

filterObj = Adapt()

# Show signal over time and signal spectrum
if 0:
    # 3.3.1
    filterObj.drawSignal(variklioSig, "variklis")
    filterObj.drawSignal(kabinosSig, "kabinos")
    filterObj.drawSignal(pilotoSig, "piloto")

    filterObj.drawSpectrum(variklioSig, "variklis")
    filterObj.drawSpectrum(kabinosSig, "kabinos")
    filterObj.drawSpectrum(pilotoSig, "piloto")

    # 3.3.2
    filterObj.saveSignalAsWav(variklioSig, f"out/variklis.wav")
    filterObj.saveSignalAsWav(kabinosSig, f"out/kabina.wav")
    filterObj.saveSignalAsWav(pilotoSig, f"out/pilotas.wav")
    print("Saved wav files")
    

# 3.4 MVK 
if 0:
    pilot_afterMVK = calculate_MVK(variklioSig, kabinosSig)

    filterObj.drawSignal(pilot_afterMVK, "piloto po MKV")
    filterObj.drawSpectrum(pilot_afterMVK, "piloto po MKV")
    filterObj.saveSignalAsWav(pilot_afterMVK, f"out/pilotasMVK.wav")
    print("Calculated and saved signal after MKV")
    

# 3.5.6 - Compare with different M
if 0:
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
        filterObj.saveSignalAsWav(sig_MVK, f"out/piloto signalas, kai M = {M}, mu = {mu_string}.wav")
    
# 3.5.7
# Filtro koeficientų skaičių 𝑀 galite keisti logaritminiu masteliu: 10,
# 20, 50, 100. Adaptacijos žingsnio 𝜇 vertę galite keisti dėsniu 0.001, 0.005, 0.01, 0.05, 0.1.
if 0:
    sig_MVK = calculate_MVK(variklioSig, kabinosSig)

    print(calculate_MSE(pilotoSig, sig_MVK))

    M = [10, 15, 20, 25, 50, 100]
    u_mu = [0.001, 0.005, 0.01, 0.05, 0.1]

    print("Calculating MSE for different M and mu combinations")
    MSE_array = np.zeros((len(M), len(u_mu)))

    for i in range(len(M)):
        for j in range(len(u_mu)):
            MSE_array[i, j] = calculate_MSE(pilotoSig, calculate_MVK(variklioSig, kabinosSig, M[i], u_mu[j]))
    print("Calculated MSE")
    print(MSE_array)

    minIdx = np.unravel_index(np.argmin(MSE_array, axis=None), MSE_array.shape)

    print(minIdx[0])
    print(minIdx[1])
    print(f"Min MSE is with: M={M[minIdx[0]]}, mu = {u_mu[minIdx[1]]}")


    plt.figure   
    for i_miu in range(len(u_mu)):
        temp_mu = u_mu[i_miu]
        plt.plot(M, MSE_array[:, i_miu], marker='D', label=f"mu = {temp_mu}")
    plt.title('MSE nuo M ir mu')
    plt.xlabel('M')
    plt.ylabel('MSE')
    plt.grid(True)
    plt.legend(loc="upper left")
    plt.savefig(filterObj.saveDir + "mse_compare",  bbox_inches='tight', pad_inches=0)
    plt.show(block=False)


# plot MSE_array over time
# Pavaizduokite adaptacijos greičio kreives (𝑀𝑆𝐸 priklausomybes nuo laiko) esant adaptacijos
# žingsnio vertėms 𝜇 = 0.001, 𝜇 = 0.01, 𝜇 = 0.1. 𝑀𝑆𝐸 galite skaičiuoti 10 ms ar ilgesniuose
# intervaluose.
if 0:
    intervalMSE_s = 20 * 10**-3
    timeMSE_ms = np.arange(0, filterObj.time_s*1000, intervalMSE_s)

    bestM = 15
    test_mu = [0.001, 0.01, 0.1]
    pilotEstimate1 = calculate_MVK(variklioSig, kabinosSig, bestM, test_mu[0])
    pilotEstimate2 = calculate_MVK(variklioSig, kabinosSig, bestM, test_mu[1])
    pilotEstimate3 = calculate_MVK(variklioSig, kabinosSig, bestM, test_mu[2])

    MSE_overTime_array = np.zeros((3, (len(timeMSE_ms))))

    start_i = 0
    end_i = int(intervalMSE_s*filterObj.Fs_Hz) # python needs int

    print("Calculating MSE over time")
    increaseInterval = end_i
    for i in range(len(timeMSE_ms)):
        MSE_overTime_array[0, i] = calculate_MSE(pilotoSig[start_i:end_i], pilotEstimate1[start_i:end_i])
        MSE_overTime_array[1, i] = calculate_MSE(pilotoSig[start_i:end_i], pilotEstimate2[start_i:end_i])
        MSE_overTime_array[2, i] = calculate_MSE(pilotoSig[start_i:end_i], pilotEstimate3[start_i:end_i])
        start_i = start_i + increaseInterval
        end_i = end_i + increaseInterval

    plt.figure
    for i in range(len(test_mu)):
        plt.plot(timeMSE_ms, MSE_overTime_array[i], label=f"mu = {test_mu[i]}")
    plt.title('MSE kitimas nuo laiko')
    plt.xlabel('t, s')
    plt.ylabel('MSE')
    plt.grid(True)
    plt.legend(loc="upper right")
    plt.savefig(filterObj.saveDir + "mse_overtime_mu",  bbox_inches='tight', pad_inches=0)
    plt.show(block=False)


# 3.6 - NMVK
# Normalizuoto mažiausių vidutinių kvadratų adaptyviojo algoritmo įgyvendinimas
if 1:
    best_mu = 0.01
    bestM = 15
    print("Comparing MSE for MVK and NMVK ")
    pilot_mvk = calculate_MVK(variklioSig, kabinosSig, bestM, best_mu)
    pilot_nmvk = calculate_NMVK(variklioSig, kabinosSig, bestM, best_mu)

    filterObj.drawSignal(pilot_mvk, "piloto naudojant MVK")
    filterObj.drawSignal(pilot_nmvk, "piloto naudojant NMVK")

    mse_mvk = calculate_MSE(pilotoSig, pilot_mvk)
    mse_mvk_mean = calculate_MSE(pilotoSig, pilot_nmvk)

    print(f"MSE of MVK: {mse_mvk}")
    print(f"MSE of MVK mean: {mse_mvk_mean}")
    
    # Save files for MATLAB:
    mvk_dic = {"pilot_mvk": pilot_mvk, "pilot_nmvk": pilot_nmvk, "data_label": "Generated signals"}
    nmvk_dic = { "label": "pilot signal nmvk"}
    scipy.io.savemat("pilot_sig.mat", mvk_dic)
    # scipy.io.savemat("pilot_sig_mvk.mat", mvk_dic)
    # scipy.io.savemat("pilot_sig_nmvk.mat", nmvk_dic)


# 4 Papildoma
if 0:
    sig_RMK = calculate_RMK(variklioSig, kabinosSig)

    filterObj.drawSignal(sig_RMK, "piloto naudojant RMK")
    print(f"RMK MSE value: {calculate_MSE(pilotoSig, sig_RMK)}")

    M = [10, 15, 20, 25, 50, 100]
    lambda_val = [1, 0.99, 0.98, 0.97, 0.96, 0.95]

    MSE_array_RMK = np.zeros((len(M), len(lambda_val)))

    print("Started calculting MSE for RMK filter")
    for i in range(len(M)):
        for j in range(len(lambda_val)):
            MSE_array_RMK[i, j] = calculate_MSE(pilotoSig, calculate_RMK(variklioSig, kabinosSig, M[i], lambda_val[j]))
    print("Calculated MSE for RMK filter")

    minIdx = np.unravel_index(np.argmin(MSE_array_RMK, axis=None), MSE_array_RMK.shape)
    print(f"Min MSE for RMK is with: M = {M[minIdx[0]]}, lambda = {lambda_val[minIdx[1]]}")
            
    plt.figure
    for i_lam in range(len(lambda_val)):
        temp_lam = lambda_val[i_lam]
        plt.plot(M, MSE_array_RMK[:, i_lam], marker='D', label=f"lambda = {temp_lam}")
    plt.title('MSE nuo M ir lamdos')
    plt.xlabel('M')
    plt.ylabel('MSE')
    plt.grid(True)
    plt.legend(loc="upper left")
    plt.savefig(filterObj.saveDir + "rmk_diff_lam", bbox_inches='tight', pad_inches=0)
    plt.show(block=False)


# compare all there metdods overtime
if 0:
    print("Comparing all 3 metdods")
    intervalMSE_s = 40 * 10**-3
    timeMSE_ms = np.arange(0, filterObj.time_s*1000, intervalMSE_s)

    bestM = 15
    bestLambda = 1
    best_mu = 0.01
    pilotEstimate1 = calculate_MVK(variklioSig, kabinosSig, bestM, best_mu)
    pilotEstimate2 = calculate_NMVK(variklioSig, kabinosSig, bestM, best_mu)
    pilotEstimate3 = calculate_RMK(variklioSig, kabinosSig, bestM, bestLambda)

    MSE_overTime_array = np.zeros((3, (len(timeMSE_ms))))

    start_i = 0
    end_i = int(intervalMSE_s*filterObj.Fs_Hz) # python needs int
    
    increaseInterval = end_i
    for i in range(len(timeMSE_ms)):
        MSE_overTime_array[0, i] = calculate_MSE(pilotoSig[start_i:end_i], pilotEstimate1[start_i:end_i])
        MSE_overTime_array[1, i] = calculate_MSE(pilotoSig[start_i:end_i], pilotEstimate2[start_i:end_i])
        MSE_overTime_array[2, i] = calculate_MSE(pilotoSig[start_i:end_i], pilotEstimate3[start_i:end_i])
        start_i = start_i + increaseInterval
        end_i = end_i + increaseInterval


    plt.figure
    plt.plot(timeMSE_ms, MSE_overTime_array[0], label="MVK")
    plt.plot(timeMSE_ms, MSE_overTime_array[1], label="NMVK")
    plt.plot(timeMSE_ms, MSE_overTime_array[2], label="RMK")
    plt.legend(loc="upper right")
    plt.title('Algoritmu palyginimas')
    plt.xlabel('t, s')
    plt.ylabel('MSE')
    plt.grid(True)
    plt.ylim(0, 0.0025)
    plt.savefig(filterObj.saveDir + "Algoritmu palyginimas",  bbox_inches='tight', pad_inches=0)
    plt.show(block=False)

 