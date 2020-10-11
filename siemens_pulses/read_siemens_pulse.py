#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:53:43 2020

@author: tangir
"""

from IPython import get_ipython
import numpy as np
import matplotlib.pyplot as plt
get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9

# %% values from PRESS paper
# 180 mao pulses duration 5.25ms, bandwidth 1.2kHz
# http://dx.doi.org/10.1002/mds.25279
pulseStr = "mao_400_4.mao_400_4.pta"
pulseDuration_ms = 5.25  # ms
n_fft = 10000

# %% Jojo's eja_svs_press protocol
pulseStr = "mao_400_4.mao_400_4.pta"
pulseDuration_ms = 7.08  # ms
n_fft = 10000

# %% Jojo's eja_svs_steam protocol
# http://dx.doi.org/10.1002/mds.25279
pulseStr = "asymm_512_6000.asymm_512_6000.pta"
pulseDuration_ms = 2.560  # ms
n_fft = 10000

# %% read/display pulse

text_file = open(pulseStr, "r")
lines = text_file.readlines()
text_file.close()

pulseShapeAmp = []
pulseShapePha = []
startReadingArray = False
for this_line in lines:
    this_line = this_line.strip()
    if(startReadingArray):
        this_line_splitted = this_line.split()
        this_amp = float(this_line_splitted[0])
        this_pha = float(this_line_splitted[1])
        pulseShapeAmp.append(this_amp)
        pulseShapePha.append(this_pha)

    if(not startReadingArray and len(this_line) == 0):
        startReadingArray = True

pulseShapeAmp = np.array(pulseShapeAmp)
pulseShapePha = np.array(pulseShapePha)
pulseShapeCmplx = pulseShapeAmp * np.exp(pulseShapePha * 1j)
n = pulseShapeCmplx.shape[0]

tvec = np.linspace(0, pulseDuration_ms / 1000.0, n)
ts = tvec[1]
fs = 1.0 / ts
fvec = np.linspace(-fs / 2.0, +fs / 2.0, n_fft)
pulseFreqResp = np.fft.fftshift(np.fft.fft(pulseShapeCmplx, n_fft))
pulseFreqRespNorm = np.abs(pulseFreqResp) / np.max(np.abs(pulseFreqResp))

# BW
fvecTop = fvec[pulseFreqRespNorm > 0.5]
bw_Hz = np.max(fvecTop) - np.min(fvecTop)
print(lines[1])
print("BW (Hz) = %fHz" % bw_Hz)
print("BWTP/R = %f" % (bw_Hz * pulseDuration_ms / 1000.0))

fig, axs = plt.subplots(2, 1)
axs[0].plot(tvec * 1000.0, np.abs(pulseShapeCmplx))
axs[0].grid('on')
axs[0].set_xlabel('time (ms)')
axs[0].set_ylabel('magnitude')

axs[1].plot(fvec, pulseFreqRespNorm)
axs[1].grid('on')
axs[1].set_xlabel('frequency (Hz)')
axs[1].set_ylabel('magnitude')



