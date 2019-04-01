'''

Analyzing PMT waveform data
20190331
Aiwu Zhang

'''

import h5py
import numpy as np
from ROOT import TCanvas, TH1F, TGraph, TH1D, TLegend, TFile, TDirectory, TTree, TStyle
import datetime
import time
from array import array
import scipy
import peakutils

if __name__ == '__main__':

    NCH = 4 # there are 4 channels
    NSamples = 10000 # each waveform has 10000 points, with a gap of 2 ns
    Nsigma = 5 # this is used to set threshold for a pulse, common for all
    # filenames to be processed
    # recommended to put three files in one process
    hdf5_file_name = ['../Data2019/D2_Ch_9_10_6_8_run5.h5', '../Data2019/D2_Ch_9_10_6_8_run6.h5', '../Data2019/D2_Ch_9_10_6_8_run7.h5']

    #binsize of charge distributions, the range is determined based on vertical range info.
    #this is supposed to be decided after data is processed and checked by user
    #(not very automatically done because there could be big difference on signal size for the four PMTs)
    BinSizeFactor = [2.5,0.5,0.5,10]

    # start a loop processing the files
    for filename in hdf5_file_name:
        # open one .h5 file the first time
        # this time I will quickly process and get two pieces of information:
        #     pulse baseline mean and width
        #     pulse amplitude and corresponding time bin information
        file = h5py.File(filename,'r')
        print filename+" is read in."

        # Final result will be saved in ROOT file
        rtfilename = filename.replace('.h5','.root')
        rootfile = TFile(rtfilename,"recreate")
        print rtfilename + " is created."

        # in the ROOT file, a few waveforms will be kept in directory "Waveforms"
        # the other histograms will be saved in "Results" directory
        waveDir = rootfile.mkdir("Waveforms")
        resultsDir = rootfile.mkdir("Results")

        # here, start getting data cards from the .h5 file
        PMTInfoDir = file['RunInfo/PMTInfo']
        ScopeInfoDir = file['RunInfo/ScopeInfo']
        WaveformDir = file['Waveforms']

        # for charge histograms, the range should be determined based on the vertical range setting during data taking
        charge_range = [0 for i in range(4)]
        for i in range(4):
            charge_range[i]=ScopeInfoDir[i]/1000.0*40.0*20.0

        # Prepare histograms to be saved in ROOT file
        hAmplitude_list = []
        hAmplitudeBin_list = []
        hWave_list = []
        hPedMean_list = []
        hPedWidth_list = []
        hFinalCharge_list = []
        hNbOfPulses_list = []
        for i in range(0,NCH,1):
            # will keep a few waveforms 
            name = "Wave_" + PMTInfoDir[i*2]# + str(i)
            hist = TH1F(name,"",10000,0,10000)
            hist.SetXTitle("Sample number")
            hist.SetYTitle("Amplitude (mV)")
            hist.SetLineColor(i+1)
            hWave_list.append(hist)
            # histogram of pulse baseline 
            name = "PedMean_" + PMTInfoDir[i*2]# + str(i)
            hist = TH1F(name,"",1000,-10,10)
            hist.SetXTitle("Amplitude (mv)")
            hist.SetYTitle("Counts")
            hist.SetLineColor(i+1)
            hPedMean_list.append(hist)
            # histogram of pulse baseline width
            name = "PedWidth_" + PMTInfoDir[i*2]# + str(i)
            hist = TH1F(name,"",1000,-10,10)
            hist.SetXTitle("Amplitude (mV)")
            hist.SetYTitle("Counts")
            hist.SetLineColor(i+1)
            hPedWidth_list.append(hist)
            # pulse charge due to fiber triggers, unit converted to fC
            name = "FinalCharge_" + PMTInfoDir[i*2]# + str(i)
            hist = TH1F(name,"",int(charge_range[i]*BinSizeFactor[i]),0,charge_range[i])
            hist.SetXTitle("Charge (fC)")
            hist.SetYTitle("Counts")
            hist.SetLineColor(i+1)
            hFinalCharge_list.append(hist)
            # pulse amplitude, or nimimum (because of negative pulse)
            name = "PulseAmplitude_" + PMTInfoDir[i*2]# + str(i)
            hist = TH1F(name,"",int(ScopeInfoDir[i]),0,ScopeInfoDir[i])
            hist.SetXTitle("Amplitude (mV)")
            hist.SetYTitle("Counts")
            hist.SetLineColor(i+1)
            hAmplitude_list.append(hist)
            # time bin where pulse amplitude is found, used to determine charge integration region
            name = "PulseAmplitudeBin_" + PMTInfoDir[i*2]# + str(i)
            hist = TH1F(name,"",10000,0,10000)
            hist.SetXTitle("Amplitude bin number")
            hist.SetYTitle("Counts")
            hist.SetLineColor(i+1)
            hAmplitudeBin_list.append(hist)
            # number of pulses found in the waveform not in the fiber trigger region
            name = "NbOfPulses_" + PMTInfoDir[i*2]# + str(i)
            hist = TH1F(name,"",20,0,20)
            hist.SetXTitle("Number of pulses from bin 1500 (3 #mus) to 10000 (20 #mus)")
            hist.SetYTitle("Counts")
            hist.SetLineColor(i+1)
            hNbOfPulses_list.append(hist)
 
        # calculate basline
        # using the first 200 sampling points
        baseline_mean = [0., 0., 0., 0.]
        baseline_width = [0., 0., 0., 0.]
        threshold = [0., 0., 0., 0.]
        eventNb = -1
        for wavekey in WaveformDir.keys():
            eventNb += 1
            thiswaveform = np.array(WaveformDir[wavekey])
            # save waveforms every 50 events
            if np.mod(int(eventNb),50)==0: #save a waveform every 100 events
                print eventNb
                for ch in range(0,NCH,1):
                    for bin in range(0,NSamples,1):
                        hWave_list[ch].SetBinContent(bin+1,thiswaveform[ch][bin]*1000.0)
                    waveDir.cd()
                    newName = "Wave_Ch_" + str(ch) + "_" + str(eventNb)
                    hWave_list[ch].SetName(newName)
                    hWave_list[ch].Write()
            # the main processing part
            for ch in range(0,NCH,1):
                baseline_mean[ch] = np.average(thiswaveform[ch][0:200])
                baseline_width[ch]= np.std(thiswaveform[ch][0:200])
                threshold[ch] = baseline_mean[ch] - Nsigma*baseline_width[ch]
                hPedMean_list[ch].Fill(baseline_mean[ch])
                hPedWidth_list[ch].Fill(baseline_width[ch])
                # look for pulse amplitude time bin in region < 1500 points
                # during data taking this should be around 600 or 1000 points
                TimeBinOfAmplitude = np.argmin(thiswaveform[ch][200:1500])
                # get charge integration about the amplitude: 20 ns before to 36 ns after the amplitude
                sumcharge = np.sum(thiswaveform[ch][TimeBinOfAmplitude+200-10:TimeBinOfAmplitude+200+19])
                fC = (sumcharge-baseline_mean[ch]*(19-10+1))*40.0 # 2 ns interval on 50 Ohm resister
                #print eventNb, ch, baseline_mean[ch], threshold[ch], thiswaveform[ch][TimeBinOfAmplitude], TimeBinOfAmplitude
                if thiswaveform[ch][TimeBinOfAmplitude+200] < threshold[ch]:
                    hAmplitude_list[ch].Fill(-1000.0*(thiswaveform[ch][TimeBinOfAmplitude+200]-baseline_mean[ch]))
                    hAmplitudeBin_list[ch].Fill(TimeBinOfAmplitude+200)
                    hFinalCharge_list[ch].Fill(-1.0*fC)
                # search for number of pulses from bin 1500 to 10000
                # that means a time window of 17 us where fiber trigger should not show but dark counts appear
                peakindex = peakutils.peak.indexes(-1.0*thiswaveform[ch][1500:NSamples],thres=abs(threshold[ch]), min_dist=40, thres_abs=True)
                #if len(peakindex)>0:
                hNbOfPulses_list[ch].Fill( len(peakindex))
                #one can print out the peak search result and check
                #print eventNb, ch, peakindex+1500

        file.close()

        #finally save all results in the ROOT file
        resultsDir.cd()
        for i in range(0,NCH,1):
            hFinalCharge_list[i].Write()
        for i in range(0,NCH,1):
            hAmplitude_list[i].Write()
        for i in range(0,NCH,1):
            hAmplitudeBin_list[i].Write()
        for i in range(0,NCH,1):
            hPedMean_list[i].Write()
        for i in range(0,NCH,1):
            hPedWidth_list[i].Write()
        for i in range(0,NCH,1):
            darkrate = hNbOfPulses_list[i].Integral(2,20)/(17.0*1.0e-6*hNbOfPulses_list[i].GetEntries())
            title = "dark rate on average " + str(darkrate) + " Hz"
            hNbOfPulses_list[i].SetTitle(title)
            hNbOfPulses_list[i].Write()
        rootfile.Close()
               

