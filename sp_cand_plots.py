import your 
from your.candidate import Candidate, crop
from your.utils.plotter import plot_h5
import numpy as np
from scipy.signal import detrend
import os
from argparse import ArgumentParser
from subprocess import call
import inspect
import matplotlib.pyplot as plt 

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

##############################
## Making and Reading Cands ##
##############################

class SP_CAND:
    def __init__(self, line, ctype):
        self.dm = None
        self.snr = None
        self.time = None
        self.samp = None
        self.wbins = None

        if ctype == 'sp':
            self.input_sp(line)

        elif ctype == 'cand':
            self.input_cand(line)
    
        else:
            print("Unrecognized data type: %s" %ctype)

    def input_sp(self, line):
        cols = line.split()
        self.dm = float(cols[0])
        self.snr = float(cols[1])
        self.time = float(cols[2])
        self.samp = int(cols[3])
        self.wbins = int(cols[4])

    def input_cand(self, line):
        cols = line.split()
        self.dm = float(cols[6])
        self.snr = float(cols[0])
        self.time = float(cols[3])
        self.samp = int(cols[1])
        self.wbins = int(cols[4])

    def __str__(self):
        str = "SP_CAND(t=%.2f s, SNR=%.2f)" %(self.time, self.snr)
        return str

    def __repr__(self):
        str = "SP_CAND(t=%.2f s, SNR=%.2f)" %(self.time, self.snr)
        return str


def cands_from_spfile(spfile, ctype='sp'):
    """
    Read in the candidates from a *singlepulse (sp) file 
    or *cand (cand) file and return an array of SP_CAND 
    class objects
    """    
    candlist = []
    with open(spfile, 'r') as fin:
        for line in fin:
            if line[0] in ["\n", "#", 'i']:
                continue
            else: pass
    
            SP = SP_CAND(line, ctype)
            candlist.append(SP)

    return np.array(candlist)




################################
## Cross Matching Cand Lists  ##
################################

def cross_check_clists(clist1, clist2):
    """
    Check for cross-matches for clist1 in clist2
    """
    # Samples 
    ts1 = np.array([ cc.samp for cc in clist1 ]) 
    ts2 = np.array([ cc.samp for cc in clist2 ]) 

    # Widths 
    w1 = np.array([ cc.wbins for cc in clist1 ])  
    w2 = np.array([ cc.wbins for cc in clist2 ])  

    # match lists 
    midx1 = []
    midx2 = []

    for ii in range(len(clist1)):
        tdiff = ts1[ii] - ts2
        jj = np.argmin( np.abs( tdiff ) )
        td_ij = np.abs(tdiff[jj])
        if td_ij <= (w1[ii] + w2[jj])//2:
            midx1.append(ii)
            midx2.append(jj)
        else:
            pass

    midx1 = np.array(midx1)
    midx2 = np.array(midx2)

    return midx1, midx2
             
        
          


#####################
## YOUR Candidates ##
#####################

def sp_to_your_cand(sp_cand, data_file, zap_chans=None, 
                    label=-1, min_samp=256, device=0):
    """
    Convert an SP_CAND object to a your Candidate obect
    
    zap_chans = list of channel numbers (data file order)
                to zero out when de-dispersing + making 
                dm time plots
    """
    cand = Candidate(fp=data_file,
                     dm=sp_cand.dm,
                     tcand=sp_cand.time,
                     width=sp_cand.wbins,
                     label=label,
                     snr=sp_cand.snr,
                     min_samp=min_samp,
                     device=device
                    )

    if zap_chans is not None:
        # Make sure chans are ints
        zz = zap_chans.astype('int')
        
        # Make kill mask array (true = flag)
        kill_mask = np.zeros(cand.nchans, dtype='bool')
        kill_mask[zz] = True
        cand.kill_mask = kill_mask
    return cand


def cand_data(cand, tstart=None, tstop=None):
    """
    Generate various data products and save your candidate 
    to h5 format 
    """
    # Get data chunk from inifle
    print("Reading candidate data")
    cand.get_chunk(tstart=tstart, tstop=tstop)

    # Make dm time data
    print("Making DM vs Time data")
    cand.dmtime()
    
    # Dedispersing
    print("Dedispersing candidate data")
    cand.dedisperse() 

    # Save candidate 
    cand.dm_opt = -1
    cand.snr_opt = -1

    return cand


def save_plot(h5file, detrend=False):
    """
    Make and save a cand plot from h5 file
    """
    plt.ioff()
    plot_h5(h5file, detrend_ft=detrend, save=True)
    plt.close()
    plt.ion()
    return 


def get_start_stop(cand, n):
    """
    Get start / stop times for data chunk

    Will take data half-length to be the larger of:

      n * max(width // 2, 1) * tsamp 

    or 

      (dispersion sweep) + width * tsamp  
    """
    w = cand.width 
    dt = cand.native_tsamp 
    tdm = cand.dispersion_delay() 
    tcand = cand.tcand

    tdur1 = n * max(w//2, 1) * dt
    tdur2 = tdm + w * dt 

    tdur = max( tdur1, tdur2 )

    tstart = tcand - tdur
    tstop  = tcand + tdur

    return tstart, tstop 


def make_h5_from_sp(sp, data_file, snr_min=7, zap_chans=None, 
                    t_dec=-1, f_dec=1, min_samp=0, n=-1, ctype='sp', out_dir='/'):
    """
    Make h5 candidates from spfile or candlist.  
    Only save candidates with snr > snr_min

    t_dec == - 1 means use width//2
    """
    # Check if candlist or sp file 
    
    # Assuming string is file name  
    if type(sp) == str:
        clist_all = cands_from_spfile(sp, ctype)

    # If list or array assume it is clist
    elif type(sp) in [np.ndarray, list]:
        clist_all = [ cc for cc in sp ] 

    # Otherwise we dont know
    else: 
        print("Unknown sp type: must be file name or array/list")
        return

    # Make new array with snr_min cut
    clist = np.array([ cc for cc in clist_all if cc.snr > snr_min ])

    # Print number of cands
    print("Processing %d cands" %len(clist))

    # now convert to your format and save
    for cc in clist:
        # convert to your format 
        yc = sp_to_your_cand(cc, data_file, zap_chans=zap_chans, 
                             min_samp=min_samp)

        # if n > 0, get tstart, tstop otherwise none
        if n > 0:
            tstart, tstop = get_start_stop(yc, n)
        else:
            tstart = tstop = None
            
        yc = cand_data(yc, tstart, tstop)

        # Decimate if desired
        t_dec_fac = t_dec
        f_dec_fac = f_dec
 
        if t_dec == -1:
            t_dec_fac = max(1, yc.width // 2)

        if t_dec_fac > 1 or f_dec_fac > 1:
            print(t_dec_fac, f_dec_fac)
            yc = decimate_cand(yc, t_dec_fac, f_dec_fac)
        else: 
            pass

        # Save datapi
        fout = yc.save_h5(out_dir=out_dir)
    return 


def decimate_cand(cand, time_fac, freq_fac):
    """
    avg in time and freq 
    """
    # Decimate time axis of tf data
    cand.decimate(key="ft", axis=0, pad=True, 
                  decimate_factor=time_fac, mode="median")

    # Decimate freq axis of tf data
    cand.decimate(key="ft", axis=1, pad=True, 
                  decimate_factor=freq_fac, mode="median")

    # Decimate dmt data
    cand.decimate(key="dmt", axis=1, pad=True, 
                  decimate_factor=time_fac, mode="median")

    return cand


def parse_input():
    """
    Use 'your' to make pulse candidates + plots
    """
    prog_desc = "Make plots of SP cand file"
    parser = ArgumentParser(description=prog_desc)

    parser.add_argument('infile', help='Input *.singlepulse file')
    parser.add_argument('filfile', help='Filterbank file')
    parser.add_argument('-s', '--snr',
                        help='SNR threshold for plotting (default = all)',
                        required=False, type=float, default=-1)
    parser.add_argument('-t', '--tdec',
                        help='Time decimation factor (default = -1, auto)',
                        required=False, type=int, default=-1)
    parser.add_argument('-f', '--fdec',
                        help='Frequency decimation factor (default = 1)',
                        required=False, type=int, default=1)
    parser.add_argument('-n', '--nmin',
                        help='Minimum output time bins (default = 32)',
                        required=False, type=int, default=32)
    parser.add_argument('-o', '--out',
                        help = 'out directory',
                        required=False, type=str, default=None)
    args = parser.parse_args()

    infile = args.infile
    filfile = args.filfile
    snr = args.snr
    tdec = args.tdec
    fdec = args.fdec
    nmin = args.nmin
    outdir = args.out

    return infile, filfile, snr, tdec, fdec, nmin, outdir

debug = 0

if __name__ == "__main__":
    if debug:
        pass
    
    else:
        infile, filfile, snr, tdec, fdec, nmin, outdir = parse_input()
        
        # Generate *h5 data snippets from cand file
        make_h5_from_sp(infile, filfile, snr_min=snr, 
                        t_dec=tdec, f_dec=fdec, n=nmin, out_dir=outdir)

        # Make plots from *h5 files using your script
        call("your_h5plotter.py -f %s/*h5 -o %s" % (outdir, outdir), shell=True)
