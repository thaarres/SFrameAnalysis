import sys, os, glob

run = True
if __name__ == "__main__":
    
    folder = sys.argv[1]
    outfolder = "HaddedOutput"
    
    if len(sys.argv) > 2:
        outfolder = sys.argv[2]

    files = dirList = glob.glob(folder+"/*root")
    counter = 0
    samples = []
    ana = []
    for f in files:
        #print f.split('.')
        sample = f.split('.')[1]
        samples.append(sample)

        
    samples = set(samples)
    print samples
    
    
    try: os.stat(outfolder) 
    except: os.mkdir(outfolder)
    
    for sample in samples:
        cmd = "hadd -f " + outfolder + "/" + sample + ".root " + folder + "/*" + sample + "*" 
        print cmd
        os.system(cmd)
        
