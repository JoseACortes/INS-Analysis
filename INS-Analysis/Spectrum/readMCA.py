def read_spectrum(file_path):
    spectrum = []
    with open(file_path, 'r') as file:    
        for line in file:
            spectrum.append(line.strip())
    return spectrum

def format_spectrum(spectrum):
    return list(map(int, spectrum[spectrum.index("2048")+1:]))

def readMCA(file_path):
    spectrum = read_spectrum(file_path)
    spectrum = format_spectrum(spectrum)
    return spectrum

    