import numpy as np
import sys

narr = np.random.randint(100, size=(4, 512))
narr = narr / 100


np.set_printoptions(threshold=sys.maxsize)

print("Real Data")
print(narr)

print()

print("FFT Data")
fnarr = np.abs(np.fft.rfft(narr))
print(fnarr)

