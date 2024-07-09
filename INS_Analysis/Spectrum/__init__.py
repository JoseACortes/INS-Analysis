from .readMCA import readMCA
from .readMCTAL import readMCTAL

filereader = {'mca': readMCA, 'mctal': readMCTAL}

def read(filename, **kwargs):
    ext = filename.split('.')[-1]
    if ext in filereader:
        return filereader[ext](filename, **kwargs)
    else:
        raise ValueError('File extension not recognized')