import h5py
import numpy as np
#import matplotlib.pyplot as plt
#import graphUtils
from ROOT import TCanvas, TH1F, TGraph, TH1D, TLegend, TFile, TDirectory, TTree
import datetime
import time
from array import array

NCH = 4 # there are 4 channels
NSamples = 10000 # each waveform has 1000 points, with a gap of 0.2 ns
hdf5_file_name = ['../Data2019/D2_Ch_9_10_6_8_run7.h5']
NSigma = 5.0 # I use 5 sigma to cut on pedestal width
rootfile = TFile("D2_Ch_9_10_6_8_run7.root","recreate")
#histograms for holding waveform pedestals
hPed_list = []
hCharge_list = []
hWave_list = []
hPedMean_list = []
hPedWidth_list = []
hFinalCharge_list = []

for i in range(0,NCH,1):
    name = "Ped_Ch_" + str(i)
    hist = TH1F(name,"",200,-5*0.001,5*0.001)
    hist.SetXTitle("Charge (mV)")
    hist.SetYTitle("Counts")
    hist.SetLineColor(i+1)
    hPed_list.append(hist)
    name = "Charge_Ch_" + str(i)
    hist = TH1F(name,"",200,0,10000)
    hist.SetXTitle("Charge (mV)")
    hist.SetYTitle("Counts")
    hist.SetLineColor(i+1)
    hCharge_list.append(hist)
    name = "Wave_Ch_" + str(i)
    hist = TH1F(name,"",1000,0,1000)
    hist.SetXTitle("Sample number")
    hist.SetYTitle("Amplitude (mV)")
    hist.SetLineColor(i+1)
    hWave_list.append(hist)
    name = "PedMean_Ch_" + str(i)
    hist = TH1F(name,"",500,-5,5)
    hist.SetXTitle("Amplitude (mv)")
    hist.SetYTitle("Counts")
    hist.SetLineColor(i+1)
    hPedMean_list.append(hist)
    name = "PedWidth_Ch_" + str(i)
    hist = TH1F(name,"",500,-5,5)
    hist.SetXTitle("Amplitude (mV)")
    hist.SetYTitle("Counts")
    hist.SetLineColor(i+1)
    hPedWidth_list.append(hist)
    name = "FinalCharge_Ch_" + str(i)
    hist = TH1F(name,"",100000,0,100000)
    hist.SetXTitle("Charge (fC)")
    hist.SetYTitle("Counts")
    hist.SetLineColor(i+1)
    hFinalCharge_list.append(hist)

waveDir = rootfile.mkdir("Waveforms")
resultsDir = rootfile.mkdir("Results")

for filename in hdf5_file_name:
    file = h5py.File(filename,'r')
    WaveformDir = file['Waveforms']
    eventNb = -1
    for wavekey in WaveformDir.keys():
        eventNb += 1
        thiswaveform = WaveformDir[wavekey]
        if np.mod(int(eventNb),100)==0:
            print eventNb
            for ch in range(0,NCH,1):
                for bin in range(0,NSamples,1):
                    hWave_list[ch].SetBinContent(bin+1,thiswaveform[ch][bin]*1000.0)
                waveDir.cd()
                newName = "Wave_Ch_" + str(ch) + "_" + str(eventNb)
                hWave_list[ch].SetName(newName)
                hWave_list[ch].Write()
        ped_mean = []
        threshold = []
        for ch in range(0,NCH,1):
            for bin in range(8000,NSamples,1):
                hPed_list[ch].Fill(thiswaveform[ch][bin])
            #
            hPed_list[ch].Fit("gaus","Q")
            mean = hPed_list[ch].GetFunction("gaus").GetParameter(1)
            sigma = hPed_list[ch].GetFunction("gaus").GetParameter(2)
            threshold = mean - NSigma*sigma
            hPedMean_list[ch].Fill(mean*1000.0)
            hPedWidth_list[ch].Fill(sigma*1000.0)
            #
            charge = 0
            maxq = 100
            fC = 0
            for bin in range(100,NSamples,1):
                if thiswaveform[ch][bin]<threshold:
                    charge += thiswaveform[ch][bin] - mean
                    fC += (thiswaveform[ch][bin]-mean)*2*1.0e6
                if maxq > thiswaveform[ch][bin]:
                    maxq = thiswaveform[ch][bin]
            if maxq < threshold:
                hCharge_list[ch].Fill(-1000.0*charge)
                hFinalCharge_list[ch].Fill(-1.0*fC)
            #
            for bin in range(0,200,1):
                hPed_list[ch].SetBinContent(bin+1,0)
#        if eventNb == 20:    break
    file.close()

resultsDir.cd()
for i in range(0,NCH,1):
    hCharge_list[i].Write()
for i in range(0,NCH,1):
    hFinalCharge_list[i].Write()
for i in range(0,NCH,1):
    hPedMean_list[i].Write()
for i in range(0,NCH,1):
    hPedWidth_list[i].Write()
rootfile.Close()

