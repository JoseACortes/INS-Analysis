try:
    from mcnptools import Mctal, MctalTally
except:
    print("warning: mcnptools is not installed, some functions will not work.")
    Mctal = None
    MctalTally = None

def readMCTAL(file, tally=8, start_time_bin=0, end_time_bin=None, nps=1e9):
    m = Mctal(file)
    tfc = MctalTally.TFC

    t8 = m.GetTally(tally)

    t8_e_bins = t8.GetEBins()
    if end_time_bin is None:
        end_time_bin = start_time_bin
    #                        f    d    u    s    m    c   e   t
    t8_evals = [[t8.GetValue(tfc, tfc, tfc, tfc, tfc, tfc, i, time) * nps for i in range(len(t8_e_bins))] for time in range(start_time_bin, end_time_bin+1)]
    
    if len(t8_evals) == 1:
        t8_evals = t8_evals[0]

    return t8_e_bins, t8_evals