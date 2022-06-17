import numpy as np
import sys
from scipy.io import wavfile
import subprocess
import os
import time
import csv

time = 2;
vec_num=8;
bin_size=vec_num*8;

print("1")
np.set_printoptions(threshold=sys.maxsize)

base_names = ["secured_conversation_", "secured_music_", "secured_cooking_audio_"]
data_directory="."
print("2")
def get_audio(directory, name):
	sr, data = wavfile.read(directory + "/" + name)
	group = data / 32767
	return np.array(group, dtype=np.float32)
print("3")
for i in range(len(base_names)):
	track1_name=base_names[i]+"48khz_track1_ds2.wav"
	track2_name=base_names[i]+"48khz_track1_ds2.wav"
	track3_name=base_names[i]+"48khz_track1_ds2.wav"
	print("4")
	track1=get_audio(data_directory, track1_name)
	track2=get_audio(data_directory, track2_name)
	track3=get_audio(data_directory, track3_name)
print("5")
print(track1)

#goal: dump data into c file, float array and time the program
#method: take timer for seconds, vec, vec_num+time_size=vec_len, bin size
#2 secs of data=lots of samples, split that into 2^8 vectors, get the closest amount of samples to 2 seconds that divides evenly into 256
#constant of 2 secs
#secure folder, pick one of the immediate ones, conversation/music/cooking audio


#use norm/other python script to extract bits
#use the output of that to pipe it into a c file
