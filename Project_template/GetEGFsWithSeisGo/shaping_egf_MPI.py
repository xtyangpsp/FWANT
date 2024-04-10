import sys,time,os
from multiprocessing import Pool
from seisgo.noise import shaping_corrdata
from seisgo import utils
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

'''
This script convolves the two sides data with a gaussian wavelet and save it.
'''
########################################
#########PARAMETER SECTION##############
########################################

# absolute path parameters
def main():
    #get arguments on the number of processors
    narg=len(sys.argv)
    if narg == 1:
        nproc=1
    else:
        nproc=int(sys.argv[1]) 

    ## Global parameters
    rootpath = "data_craton"                                 # root path for this data processin
    INDIR  = os.path.join(rootpath,'PAIRS_TWOSIDES')                          # dir where stacked data is going to
    
    wavelet="gaussian"
    wave_width=2
    wave_shift=3*wave_width 
    stack=True 
    trim_end=True
    stack_method="robust"
    comp="ZZ"
    output_format="asdf"  #"sac" or "asdf"

    OUTDIR_ROOT    = INDIR+"_"+wavelet+"_a"+str(wave_width)+"t"+str(wave_shift)
    if not os.path.isdir(OUTDIR_ROOT):os.makedirs(OUTDIR_ROOT)
    #######################################
    ###########PROCESSING SECTION##########
    #######################################

    #loop through resources, for each source MPI through station pairs.
    sources_temp=utils.get_filelist(INDIR)
    #exclude non-directory item in the list
    sources=[]
    for src in sources_temp:
        if os.path.isdir(src): sources.append(src)
    if nproc >=2:
        p=Pool(int(nproc))
    tt0=time.time()
    for src in sources:
        # cross-correlation files
        ccfiles = utils.get_filelist(src,"h5")
        print("assembled %d files"%(len(ccfiles)))
        outdir=os.path.join(OUTDIR_ROOT,os.path.split(src)[1])

        #loop for each station pair
        print("working on all pairs with %d processors."%(nproc))
        if nproc < 2:
            for j in range(len(ccfiles)):
                shaping_corrdata(ccfiles[j],wavelet,wave_width,wave_shift,trim_end,outdir,
                comp,stack,stack_method,output_format=output_format)
        else: 
            p.starmap(shaping_corrdata,[(ccfile,wavelet,wave_width,wave_shift,trim_end,outdir,\
                comp,stack,stack_method,output_format) for ccfile in ccfiles])

    print('it took %6.2fs %s' % (time.time()-tt0,src))
        # 
    if nproc >=2:
        p.close()
if __name__ == "__main__":
    main()
