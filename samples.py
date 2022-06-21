import numpy as np
import wave
import contextlib
import sys
from scipy.io import wavfile
import subprocess
import os
import time
import csv
import wave

#512x256 = 1,310.72=131,072/48=2.744 seconds of data to take, 
vec_num=256;
vec_len=512;

np.set_printoptions(threshold=sys.maxsize)

base_names = ["secured_conversation_", "secured_music_", "secured_cooking_audio_"]
data_directory="."

def get_audio(directory, name):
	sr, data = wavfile.read(directory + "/" + name)
	group = data / 32767
	return np.array(group, dtype=np.float32)

for i in range(len(base_names)):
	track1_name=base_names[i]+"48khz_track1_ds2.wav"
	track2_name=base_names[i]+"48khz_track1_ds2.wav"
	track3_name=base_names[i]+"48khz_track1_ds2.wav"
	
	track1=get_audio(data_directory, track1_name)
	track2=get_audio(data_directory, track2_name)
	track3=get_audio(data_directory, track3_name)
	
	group=np.array(track1+track2, dtype=np.float32)

#for elem in range(0,len(track1),131072):
converted_list = [str(element) for element in track1]
joined_string=", ".join(converted_list)
	

file_name='data.c'
f=open(file_name, 'w+')


elem=131073
it=101
#change it everytime u need to run a new set of data
#move data.c to embedded 

with open(file_name, 'w+') as f:
	f.write('#include "system.h"\n#include "arm_math.h"\n#define ARMCM4_FP')
	f.write("\nconst float32_t input_data[] = {")
	f.write(joined_string[it*elem:(it+1)*elem])
	f.write("}; \n")

#txt file with array, float32 array, write to c file
#print(joined_string[0:131072])

#print(track1[0:131072, sep = ", ")

#goal: dump data into c file, float array and time the program
#method: take timer for seconds, vec, vec_num+time_size=vec_len, bin size
#2 secs of data=lots of samples, split that into 2^8 vectors, get the closest amount of samples to 2 seconds that divides evenly into 256
#constant of 2 secs
#secure folder, pick one of the immediate ones, conversation/music/cooking audio


#use norm/other python script to extract bits
#use the output of that to pipe it into a c file
