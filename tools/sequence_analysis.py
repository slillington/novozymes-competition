"""
Amino acid attributes based on tables in amino_acid_scales.pptx
Sally Jiao
"""
def get_aa_size(aa_code):
    size = {'A':-0.733, 'C':-0.862, 'D':-3.656, 'E':1.477, 'F':1.891, 'G':1.330, 'H':-1.673, 'I':2.131, 'K':0.533, 'L':-1.505, 'M':2.219, 'N':1.299, 'P':-1.628, 'Q':-3.005, 'R':1.502, 'S':-4.760, 'T':2.213, 'V':-0.544, 'W':0.672, 'Y':3.097}
    return size[aa_code]

def get_aa_volume(aa_code):
    volume = {'A':92, 'R':225, 'N':135, 'D':125, 'C':106, 'Q':161, 'E':155, 'G':66, 'H':167, 'I':169, 'L':168, 'K':171, 'M':171, 'F':203, 'P':129, 'S':99, 'T':122, 'W':240, 'Y':203, 'V':142}
    return volume[aa_code]

def get_aa_alphapropensity(aa_code):
    alpha = {'A':825, 'R':297, 'N':-379, 'D':-304, 'C':-42, 'Q':413, 'E':642, 'G':-1050, 'H':-118, 'I':144, 'L':766, 'K':242, 'M':265, 'F':91, 'P':-895, 'S':-373, 'T':-391, 'W':87, 'Y':-53, 'V':-172}
    return alpha[aa_code]

def get_aa_betapropensity(aa_code):
    beta = {'A':-357, 'L':174, 'E':-392, 'Q':-247, 'R':-107, 'M':-19, 'K':-264, 'V':1280, 'I':945, 'Y':470, 'F':458, 'T':281, 'W':157, 'G':-492, 'N':-450, 'P':-664, 'D':-623, 'S':-203, 'C':94, 'H':2}
    return beta[aa_code]

def get_aa_coilpropensity(aa_code):
    coil = {'A':-248, 'L':-413, 'E':-358, 'Q':-159, 'R':-106, 'M':4, 'K':-101, 'V':-281, 'I':-326, 'Y':-184, 'F':-163, 'T':287, 'W':-89, 'G':73, 'N':217, 'P':1160, 'D':400, 'S':334, 'C':70, 'H':34}
    return coil[aa_code]
    
def get_aa_polarity(aa_code):
    polarity = {'I':-1.13, 'L':-1.18, 'V':-1.27, 'A':0.10, 'G':0.33, 'P':0.73, 'F':-2.12, 'M':-1.59, 'W':-0.51, 'Y':-0.21, 'H':-0.50, 'T':0.07, 'S':0.52, 'N':0.48, 'Q':0.95, 'D':0.78, 'E':0.83, 'K':1.40, 'R':1.91, 'C':-1.42}
    return polarity[aa_code]

def get_aa_charge(aa_code):
    charge = {'R':1, 'H':0, 'K':1, 'D':-1, 'E':-1, 'S':0, 'T':0, 'N':0, 'Q':0, 'C':0, 'G':0, 'P':0, 'A':0, 'V':0, 'I':0, 'L':0, 'M':0, 'F':0, 'Y':0, 'W':0}
    return charge[aa_code]
