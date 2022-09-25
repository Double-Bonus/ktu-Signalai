# student number = 10
# accord = Dm

# f1 = 0 # is not used

import numpy as np
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
for i, note in enumerate(notes):
    if debug:
        print(i)
        print(note)
    N.append(round(Fd / note))    

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


    X_final.append(np.concatenate([x_random, x_zeros]))
    if debug:
        print(X_final[i])



print(np.shape(X_final))
# static assert ar len == (fd*ts=44100*3=132300)
