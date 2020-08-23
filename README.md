# Airgun_Process
Batch Process Procedures of Airgun Excitation Datas for Kiwi Tool 

Usage:
  1. PreProcess.py: Process the seismic data recorded from the airgun excitations in the time period you choose .
    python PreProcess.py [year] [min_day] [max_day]
  2. Stack.py: Stack the selected seismic data according to the SNR level (e.g. every 20 stacked records) and deconvolute the stacked datas in every day from the specific year.
    python Stack.py 2015 20
