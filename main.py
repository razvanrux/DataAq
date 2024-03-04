#import numpy as np

#from scipy.io import wavfile
#samplerate, data = wavfile.read('Wav20.wav')
#wavFileInfo = open("wafeInfo.txt", "a")
#wavFileInfo.write(str(samplerate)+'\n')
#wavFileInfo.write(str(data.size)+'\n')
#wavFileInfo.close()

#print(samplerate)
#print(data.size)
#print (data)

#np.savetxt("waveData.txt", data, fmt="%2.0f")

import numpy as np
from scipy.signal import hilbert;
signal=[];
file=open("waveData.txt","r");
for line in file:
    signal.append(int(line));

analytic_signal=hilbert(signal);
amplitude_envelope=np.abs(analytic_signal);
file.close();
file=open("anvelopa.txt","w");
file.write("");
file.close();
file=open("anvelopa.txt","a");
for i in amplitude_envelope:
    file.write(str(i)+"\n");

