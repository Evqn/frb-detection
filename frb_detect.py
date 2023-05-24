import numpy as np
import mmap
import struct
import time
import os
import glob
from multiprocessing.pool import ThreadPool as Pool
import multiprocessing
from numba import jit
from subprocess import call
import filterbank as fb_old
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

def file_conversion(top_dir, chunk_dir_index, dada1, dada2):
    """
    Converts telescope data to filterbank file
    """
    tstart = time.time()
    # make directory organized by chunk index
    chunk_dir_name = '%s/chunk%d' % (top_dir, chunk_dir_index)
    os.mkdir(chunk_dir_name)
    
    print('Creating new directory: %s' % chunk_dir_name)
    # symlink both dada file

    print('Creating symlinks')
    os.symlink(dada1, chunk_dir_name + '/' + os.path.basename(dada1))
    os.symlink(dada2, chunk_dir_name + '/' + os.path.basename(dada2))

    # convert both files
    convert_cmd = 'mkfb -b -p /data/pulsars -s /data/scan.table -t xka- /data/conf_xka %s %s -n 1' % (chunk_dir_name, chunk_dir_name) 
    print(convert_cmd)
    call(convert_cmd, shell=True)
    return time.time() - tstart

def get_data(fmm, fb, nspec=-1):
    """
    read in data from memory mapped filterbank file

    if nspec = -1, read all data 

    else: read in nspec timesameples 
          (so, nchan * nspec data points)

    We are assuming that data is in floats!
    (May want to change this later)

    Assumes that we are at the right spot in the file. 
    If not sure, try a fmm.seek( nbytes ) to get where 
    you want
    """
    # Figure out how many spectra we are reading in
    if nspec <= 0:
        nspec = fb.nspec
    else: pass

    nchan = fb.nchan
    count = nspec * nchan 
    dsize = 4 # bytes per float

    # fmm.read operates on bytes
    # frombuffer allows us to convert a byte array to floats 
    dat = np.frombuffer( fmm.read( count * dsize ), 
                         dtype='float32', count=count)

    #print(dat.flags)

    # Now let's reshape our data to (nsamps, nchans)
    # The mutliply by 1.0 is a very dumb hack to make 
    # sure dat is a copy (and therefore writable?). 
    # Not really sure why it is read-only at first
    dat = np.reshape(dat, (-1, nchan)) * 1.0

    #print(dat.flags))
    return dat


def mmap_bandpass(fmm, fb):
    """
    fmm = memory mapped file object
    fb  = filterbank header class object
    """
    tstart = time.time()
    # Print what we're doing
    print("Bandpass")

    # Skip past the header to get to start of the data 
    # then read in spectral data from file, unpack from 
    # bytes to floats, and reshape 
    fmm.seek(fb.header_size)
    data = get_data(fmm, fb)

    # data is in shape (nsamps, nchans) where 
    # nsamps is number of time samples and nchans 
    # is number of frequency channels

    # Bandpass is then just the time average of 
    # the array.  However, we will first calculate  
    # the sum along this axis
    bpass = np.sum(data, axis=0)

    # Now calculate number of nonzero samples 
    # per channel.
    n_pos = np.sum( data > 0, axis=0).astype('float32')

    # Now we can get a mean that is not affected by 
    # time samples that may have been zeroed out in 
    # RFI flagging or some other reason
    xx = np.where( n_pos == 0)[0]
    if len(xx):
        n_pos[xx] = 1
    bpass /= n_pos  

    # Check if any bpass values are 0, so we 
    # don't divide by zero
    xx = np.where( bpass == 0 )[0]
    if len(xx):
        bpass[xx] = 1
    else: pass

    # Now we can finally do the bandpass correction 
    # and divide data by bpass
    data /= bpass 

    # And thats it so we can update the mem mapped file

    # Go back to start of data then convert data array 
    # to bytes and write to memory mapped file
    fmm.seek(fb.header_size)
    fmm.write( data.tobytes() )

    tstop = time.time()
    dt = tstop - tstart

    return bpass, dt


def mmap_rfi_clip(fmm, fb, nsig=5):
    """
    fmm = memory mapped file object
    fb  = filterbank header class object
    """
    tstart = time.time()
    # Print what we're doing
    print("RFI Clipping")

    # Skip past the header to get to start of the data 
    # then read in spectral data from file, unpack from 
    # bytes to floats, and reshape 
    fmm.seek(fb.header_size)
    data = get_data(fmm, fb)

    # Calculate median and std dev for each channel
    # EDIT: Though not ideal, using the mean instead 
    # of median improves speed from ~25 sec to ~9 sec.
    #meds = np.median(data, axis=0)
    meds = np.mean(data, axis=0)
    sigs = np.std(data, axis=0)

    # Where does data exceed meds + nsig * sigs
    xx = np.where( data > meds + nsig * sigs )

    # How much was flagged
    print(len(xx[0]))
    flag_frac = len(xx[0]) / data.size
    print("Flagging %.2f%% of data" %(flag_frac * 100))

    # Mask them by setting to zero
    data[xx] = 0

    # And thats it so we can update the mem mapped file

    # Go back to start of data then convert data array 
    # to bytes and write to memory mapped file
    fmm.seek(fb.header_size)
    fmm.write( data.tobytes() )

    tstop = time.time()
    dt = tstop - tstart

    return dt


def mmap_rfi_narrow(fmm, fb, nsig=5):
    """
    fmm = memory mapped file object
    fb  = filterbank header class object
    """
    tstart = time.time()
    # Print what we're doing
    print("RFI Narrow band")

    # Skip past the header to get to start of the data 
    # then read in spectral data from file, unpack from 
    # bytes to floats, and reshape 
    fmm.seek(fb.header_size)
    data = get_data(fmm, fb)

    # Calculate the std dev for each channel
    sigs = np.std(data, axis=0)

    # Calc stats of std devs
    s_med = np.median(sigs)
    s_sig = np.std(sigs)

    # Find anomalous channels
    #xx = np.where( sigs > s_med + nsig * s_sig )[0]
    xx = np.where( np.abs(sigs - s_med) > nsig * s_sig )[0]

    # How much was flagged
    print("Flagging %d channels" %len(xx))

    # Mask them by setting to zero
    if len(xx):
        data[:, xx] = 0 
    else: pass

    # And thats it so we can update the mem mapped file

    # Go back to start of data then convert data array 
    # to bytes and write to memory mapped file
    fmm.seek(fb.header_size)
    fmm.write( data.tobytes() )

    tstop = time.time()
    dt = tstop - tstart

    return dt


def mmap_meansub_norm(fmm, fb):
    """
    fmm = memory mapped file object
    fb  = filterbank header class object
    """
    tstart = time.time()
    # Print what we're doing
    print("Subtract mean and normalize to unit variance")

    # Skip past the header to get to start of the data 
    # then read in spectral data from file, unpack from 
    # bytes to floats, and reshape 
    fmm.seek(fb.header_size)
    data = get_data(fmm, fb)

    # Calculate mean and the std dev for each channel
    avgs = np.mean(data, axis=0)
    sigs = np.std(data, axis=0)

    # If sigs is 0, set to 1 so we dont get divide errors
    sigs[ sigs == 0 ] = 1

    mdat = np.zeros(data.shape, dtype='float32') + avgs 
    mdat[ data == 0 ] = 0 

    # subtract mean and divide by standard deviation
    data = (data - mdat) / sigs

    # And thats it so we can update the mem mapped file

    # Go back to start of data then convert data array 
    # to bytes and write to memory mapped file
    fmm.seek(fb.header_size)
    fmm.write( data.tobytes() )

    tstop = time.time()
    dt = tstop - tstart

    return dt


def write_outfile(outfile, fmm, fb):
    """
    Write data to output file outfile
    """
    tstart = time.time()
    fmm.seek(0)
    with open(outfile, 'wb') as fout:
        fout.write( fmm.read() )

    tstop = time.time()
    dt = tstop - tstart

    return dt

def sumpols(fmm1, fb1, fmm2, fb2):
    """
    Sums two polarizations and writes to first memory mapped file
    """
    tstart = time.time()
    # Skip past headers and get data for both files
    fmm1.seek(fb1.header_size)
    data1 = get_data(fmm1, fb1)
    fmm2.seek(fb2.header_size)
    data2 = get_data(fmm2, fb2)

    data1 = data1 + data2

    fmm1.seek(fb1.header_size)
    fmm1.write(data1.tobytes())

    fmm2.close()

    tstop = time.time()
    dt = tstop - tstart

    return dt


def dedispersion(infile, dm):
    """
    Runs dedispersion
    """
    tstart = time.time()
    
    ts = infile.strip(".fil")

    # call prepfil
    prep_cmd = "prepdata " +\
               "-filterbank "+\
               "-nobary " +\
               "-dm " +\
               "%s " %dm +\
               "-o " +\
               "%s " %ts +\
               "%s " %infile
    print(prep_cmd)
    call(prep_cmd, shell=True)
    
    tstop = time.time()
    tdur = tstop - tstart
    
    de_out = ts + ".dat"
    return de_out, tdur


def single_pulse_search(infile, maxwidth=1.0, snr_min=8):
    """
    Runs single pulse search program to detect pulses
    
    """
    tstart = time.time()

    m = maxwidth
    t = snr_min

    #call single_pulse_search
    pulse_cmd = "single_pulse_search.py " +\
               "-m %.2f " %m +\
               "-t %.1f " %t +\
               "%s "%infile
    print(pulse_cmd)
    call(pulse_cmd, shell=True)

    tstop = time.time()
    tdur = tstop - tstart

    return tdur


def mmap_proc(inbase):
    """
    Test example to memory map a file 
    and do bandpass
    """ 
    infile = inbase
    outfile = inbase + ".corr"

    tstart = time.time()
    # Make filterbank class object to access the 
    # header and various metadata 
    fb = fb_old.FilterbankFile(infile)

    # Open filterbank file (in read/write mode) 
    # and memory map it
    fin = open(infile, "r+b")
    fmm = mmap.mmap(fin.fileno(), 0)
    #fmm.flush()

    # Do the clipping 
    print("Clipping...")
    dt_rfi1 = mmap_rfi_clip(fmm, fb, nsig=5)

    # Calculate and apply the bandpass
    print("Bandpass correction...")
    bpass, dt_bp = mmap_bandpass(fmm, fb)  
    
    # Do the narrow band channel masking
    print("Bad channel masking...")
    dt_rfi2 = mmap_rfi_narrow(fmm, fb, nsig=3)

    # Mean subtract and normalize
    print("Subtracting mean and normalizing..")
    dt_norm = mmap_meansub_norm(fmm, fb)

    # Write memmap file back to disk
    # print("Writing output file...")
    # dt_tw = write_outfile(outfile, fmm, fb)
    
    # Dedispersion
    # print("Dedispersion...")
    # de_out, dt_dd = dedispersion(outfile)

    # Single pulse search
    # print("Single Pulse Search...")
    # dt_sps = single_pulse_search(de_out)
     
    # Close mem map
    # fmm.close()

    # Close open file
    fin.close()

    tstop = time.time()
    dt = tstop - tstart

    print("RFI clipping    = %.1f sec" %dt_rfi1)
    print("BP time         = %.1f sec" %dt_bp)
    print("RFI bad chans   = %.1f sec" %dt_rfi2)
    print("Mean sub + norm = %.1f sec" %dt_norm)
    # print("Write file      = %.1f sec" %dt_tw)
    # print("Dedispersion    = %.1f sec" %dt_dd)
    # print("Pulse Search    = %.1f sec" %dt_sps)
    print("Total time      = %.1f sec" %dt)

    return fmm, fb


def multiprocess(file_bases):
    """
    Runs pipeline up to mean sub + norm on two multiple files at once
    """

    tstart = time.time()

    """# stores return values
    manager = mp.Manager()
    return_dict = manager.dict()

    #holds proceses to be joined
    processes = []

    for i in range(len(file_bases)):
        p = mp.Process(target = mmap_proc, args=(file_bases[i], i, return_dict))
        p.start()
        processes.append(p)
    
    for process in processes:
        process.join()"""

    p = Pool(processes = 2)
    data = p.map(mmap_proc, file_bases)
    p.close()

    tstop = time.time()
    dt = tstop - tstart

    return data, dt


def post_mp(data, outfile, dm, maxwidth=1.0, snr_min=8.0):
    """
    Handles all operations after multiprocessing
    """
    tstart = time.time()

    fmm1 = data[0][0]
    fb1 = data[0][1]
    fmm2 = data[1][0]
    fb2 = data[1][1]

    # sum polarizations
    dt_sp = sumpols(fmm1, fb1, fmm2, fb2)

    # write to outfile
    dt_tw = write_outfile(outfile, fmm1, fb1)

    # Runs dedispersion
    de_out, dt_dd = dedispersion(outfile, dm)

    # Runs single pulse search
    dt_sps = single_pulse_search(de_out, maxwidth, snr_min)
    
    print("Polarization Sum     = %.1f sec" % dt_sp)
    print("Outfile Write        = %.1f sec" % dt_tw)
    print("Dedispersion         = %.1f sec" % dt_dd)
    print("Pulse Search         = %.1f sec" % dt_sps)

    dt = time.time() - tstart
    return dt
    tstart = time.time()

def candidate_plot(chunk_dir, snr, num_chans):
    """
    Creates plots of candidates
    """
    print("Plotting candidates")

    plot_cmd = 'python /src/sp_cand_plots.py -s %d -f %d /data/chunk%d/chunk%d.singlepulse /data/chunk%d/chunk%d.fil -o /data/chunk%d/' % (snr, num_chans, chunk_dir, chunk_dir, chunk_dir, chunk_dir, chunk_dir)
    print(plot_cmd)
    call(plot_cmd, shell=True)

    return

def pipeline(chunk_dir_index, dada1, dada2):
    """
    Runs pipeline
    """
    tstart = time.time()
    top_dir = '/data'
    
    #file_bases = ["%s/m81_lcp" %top_dir, 
                  #"%s/m81_rcp" %top_dir]
    #file_bases = glob.glob('/data/*.X-LCP.fil', recursive=True) + glob.glob('/data/*.X-RCP.fil', recursive=True)
    #print(file_bases)
    outfile = "%s/chunk%d/chunk%d.fil" % (top_dir, chunk_dir_index, chunk_dir_index)

    dm        = 87.77
    maxwidth  = 1.0 
    snr_min   = 8.0
    num_chans = 10.0
    
    dt_fc = file_conversion(top_dir, chunk_dir_index, dada1, dada2)
    
    # find correct filterbank files 

    file_bases = glob.glob('/data/chunk%d/*/*X-LCP.fil' % chunk_dir_index, recursive=True) + glob.glob('/data/chunk%d/*/*X-RCP.fil' % chunk_dir_index, recursive=True)
    print(file_bases)
    
    data, dt_mp = multiprocess(file_bases)

    dt_pmp = post_mp(data, outfile, dm, maxwidth, snr_min)

    candidate_plot(chunk_dir_index, snr_min, num_chans) 

    dt = time.time() - tstart
    print("Total Time Summary")
    print("--------------------------------------")
    print("Files Conversion     = %.1f sec" % dt_fc)
    print("Multiprocessing      = %.1f sec" % dt_mp)
    print("Post-Multiprocessing = %.1f sec" % dt_pmp)
    print("Total Time           = %.1f sec" % dt)
    print("--------------------------------------")



class FileMonitor(FileSystemEventHandler):

    FILE_THRESHOLD = 1280004096
    CHUNK_SIZE = 2
    curr_file_index = 0
    process_num = 0
    flag = True


    def on_modified(self, event):
        if os.path.getsize(event.src_path) >= self.FILE_THRESHOLD and '.dada' in event.src_path:
            #print('File is done writing: %s'%event.src_path)
            # search for .dada files
            #print('File size: %d' %os.path.getsize(event.src_path))
            dada_files = sorted(glob.glob('/data/*.dada'))
            #print(dada_files)

            # make sure there are two to processes as a chunk
            while(self.curr_file_index <= len(dada_files) - self.CHUNK_SIZE):
           
                dada1 = dada_files[self.curr_file_index]
                dada2 = dada_files[self.curr_file_index+1]

                print('Multiprocessing %s and %s' % (dada1, dada2))
                # start multiprocessing
                self.process_num+=1
                p = multiprocessing.Process(target=pipeline, args=(self.process_num, dada1, dada2))
                p.start()

                self.curr_file_index += self.CHUNK_SIZE


if __name__ == "__main__":
    
    top_dir = '/data'
    interval = (5)
    #pipeline()
    event_handler = FileMonitor()
    observer = Observer()
    observer.schedule(event_handler, path=top_dir, recursive=True)
    print("Monitoring for finished file...")
    observer.start()
    try:
        while(True):
            time.sleep(interval)
    
    except KeyboardInterrupt:
        observer.stop()
        observer.join()
        
    
