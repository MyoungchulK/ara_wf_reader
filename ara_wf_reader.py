import numpy as np
import os, sys
import h5py
from matplotlib import pyplot as plt

# cern root lib
import ROOT


#input information 
#change this information by your desired run and event
Data = '/data/exp/ARA/2013/filtered/full2013Data/ARA02/root/run1449/event1449.root' # ara raw data location
Ped = '/data/exp/ARA/2013/calibration/pedestals/ARA02/pedestalValues.run001392.dat' # corresponding pedestal location. always run number should be close to 'Data' run number. And it should be smaller than 'Data' 
Station = 2 # station number
Run = 1449 # run number
Evt = 452 # event number
Output = '/data/user/mkim/OMF_sky/ARA02/' # location that you want to save the files


# import cern root and ara root lib from cvmfs
ROOT.gSystem.Load(os.environ.get('ARA_UTIL_INSTALL_DIR')+"/lib/libAraEvent.so")


### load raw data and process to general quality cut by araroot ###
# open a data file
file = ROOT.TFile.Open(Data)

# load in the event free for this file
eventTree = file.Get("eventTree")

# set the tree address to access our raw data type
rawEvent = ROOT.RawAtriStationEvent()
eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))

# get the number of entries in this file
num_events = eventTree.GetEntries()
print('total events:', num_events)

# open a pedestal file
calibrator = ROOT.AraEventCalibrator.Instance()
calibrator.setAtriPedFile(Ped, Station)

# open general quilty cut
qual = ROOT.AraQualCuts.Instance()
### load raw data and process to general quality cut by araroot ###


### convert ADC count to calibrated WF ###
# get the desire event
eventTree.GetEntry(Evt)

print('Selected event:',Evt)

# make a useful event -> calibration process
usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kLatestCalib)

# folks in your lab might want to see what is the trigger of this event:)
if rawEvent.isSoftwareTrigger() == 1:
    print('This is Software triggered event.')
elif rawEvent.isCalpulserEvent() == 1:
    print('This is Calpulser triggered event.')
elif rawEvent.isSoftwareTrigger() == 0 and rawEvent.isCalpulserEvent() == 0:
    print('This is RF triggered event.')

# make a space(list) for calibrated WF
raw_graph = []

# make a space(list) for snr
snr = []

# extracting time and volt from every antenna
for c in range(16):

    # get Tgraph(root format) for each antenna
    graph = usefulEvent.getGraphFromRFChan(c)

    # get value from Tgraph
    x_buff = graph.GetX()
    y_buff = graph.GetY()

    # into numpy array
    raw_time = np.frombuffer(x_buff,dtype=float,count=-1) # It is ns(nanosecond)
    raw_volt = np.frombuffer(y_buff,dtype=float,count=-1) # It is mV

    # put the each antennas time and volt into 'raw_graph'
    raw_graph.append(np.stack([raw_time, raw_volt],axis=-1))

    # put the each antennas snr into 'snr'
    rms = np.nanstd(raw_volt)
    peak = np.nanmax(np.abs(raw_volt))
    snr.append(peak/rms)

print('Calibration is done!')
### convert ADC count to calibrated WF ###


### save calibrated WF into h5df ###
if not os.path.exists(Output): #check whether there is a directory for save the file or not
    os.makedirs(Output) #if not, create the directory 
os.chdir(Output) #go to the directory

# create output file
h5_file_name='RawWF_ARA'+str(Station)+'_Run'+str(Run)+'_Evt'+str(Evt)+'.h5'
hf = h5py.File(h5_file_name, 'w')
    
# saving each wf seperatly
for c in range(16):
    hf.create_dataset('Ch'+str(c), data=raw_graph[c], compression="gzip", compression_opts=9)

# save snr
hf.create_dataset('SNR', data=np.asarray(snr), compression="gzip", compression_opts=9)

# close output file
hf.close()

print('output is',Output+h5_file_name)
### save calibrates WF into h5df ###


### loading wf from output and make plot ###
os.chdir(Output) #go to the directory

# load file
file = h5py.File(h5_file_name, 'r')

wf_title = 'RawWF, ARA'+str(Station)+', Run'+str(Run)+', Event'+str(Evt)
wf_file = Output + 'RawWF_ARA'+str(Station)+'_Run'+str(Run)+'_Event'+str(Evt)+'.png'

fig = plt.figure(figsize=(24, 18)) # figure size

ax = fig.add_subplot(111) # want to have a one big XY label for all 16 plots
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_xlabel(r'Time [ $ns$ ]', labelpad=30,fontsize=40)
ax.set_ylabel(r'Amplitude [ $mV$ ]',labelpad=50, fontsize=40)

plt.title(wf_title, y=1.04,fontsize=35)

# calling snr
snr = file['SNR'][:]

for b in range(16):
    # calling WF from output
    wf = file['Ch'+str(b)][:]

    # make individual wf by for loop
    ax = fig.add_subplot(4,4,b+1)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=20)
    #ax.set_xlim(-200,800)
    #ax.set_ylim(-0.8,0.8)
    ax.grid(linestyle=':')
    ax.set_title('Ch:'+str(b),fontsize=25)

    ax.plot(wf[:,0],wf[:,1],'-',lw=2,color='red',alpha=0.7,label='SNR:'+str(np.round(snr[b],1)))

    plt.legend(loc='best',numpoints = 1 ,fontsize=15)

plt.tight_layout()

# saving png into output path
fig.savefig(wf_file,bbox_inches='tight')
plt.close()
print('output is',wf_file)
### loading wf from output and make plot ###

print('Done!!')



