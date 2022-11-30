"""
Some amino acid rankings roughly correlated with various amino acid attributes
From Tables S5 and S6 of https://doi.org/10.1016/j.csbj.2020.02.021
Sally Jiao
"""
def get_aa_polarity(aa_code):
    polarity = {'A':-0.591, 'C':-1.343, 'D':1.050, 'E':1.357, 'F':-1.006, 'G':-0.384, 'H':0.336, 'I':-1.239, 'K':1.831, 'L':-1.019, 'M':-0.663, 'N':0.945, 'P':0.189, 'Q':0.931, 'R':1.538, 'S':-0.228, 'T':-0.032, 'V':-1.337, 'W':-0.595, 'Y':0.260}
    return polarity[aa_code]

def get_aa_secondarystructure(aa_code):
    ss = {'A':-1.302, 'C':0.465, 'D':0.302, 'E':-1.453, 'F':-0.59, 'G':1.652, 'H':-0.417, 'I':-0.547, 'K':-0.561, 'L':-0.987, 'M':-1.524, 'N':0.828, 'P':2.081, 'Q':-0.179, 'R':-0.055, 'S':1.399, 'T':0.326, 'V':-0.279, 'W':0.009, 'Y':0.830}
    return ss[aa_code]

def get_aa_volume(aa_code):
    vol = {'A':-0.733, 'C':-0.862, 'D':-3.656, 'E':1.477, 'F':1.891, 'G':1.330, 'H':-1.673, 'I':2.131, 'K':0.533, 'L':-1.505, 'M':2.219, 'N':1.299, 'P':-1.628, 'Q':-3.005, 'R':1.502, 'S':-4.760, 'T':2.213, 'V':-0.544, 'W':0.672, 'Y':3.097}
    return vol[aa_code]

def get_aa_codondiversity(aa_code):
    cd = {'A':1.570, 'C':-1.020, 'D':-0.259, 'E':0.113, 'F':-0.397, 'G':1.045, 'H':-1.474, 'I':0.393, 'K':-0.277, 'L':1.266, 'M':-1.005, 'N':-0.169, 'P':0.421, 'Q':-0.503, 'R':0.440, 'S':0.670, 'T':0.908, 'V':1.242, 'W':-2.128, 'Y':-0.838}
    return cd[aa_code]

def get_aa_charge(aa_code):
    charge = {'A':-0.146, 'C':-0.255, 'D':-3.242, 'E':-0.837, 'F':0.412, 'G':2.064, 'H':-0.078, 'I':0.816, 'K':1.648, 'L':-0.912, 'M':1.212, 'N':0.933, 'P':-1.392, 'Q':-1.853, 'R':2.897, 'S':-2.647, 'T':1.313, 'V':-1.262, 'W':-0.184, 'Y':1.512}
    return charge[aa_code]

def get_aa_hydrophobicity(aa_code):
    hydrophobicity = {'A':0.008, 'R':0.171, 'N':0.255, 'D':0.303, 'C':-0.132, 'Q':0.149, 'E':0.221, 'G':0.218, 'H':0.023, 'I':-0.353, 'L':-0.267, 'K':0.243, 'M':-0.239, 'F':-0.329, 'P':0.173, 'S':0.199, 'T':0.068, 'W':-0.296, 'Y':-0.141, 'V':-0.274}
    return hydrophobicity[aa_code]

def get_aa_length(aa_code):
    length = {'A':0.134, 'R':-0.361, 'N':0.038, 'D':-0.057, 'C':0.174, 'Q':-0.184, 'E':-0.280, 'G':0.562, 'H':-0.177, 'I':0.071, 'L':0.018, 'K':-0.339, 'M':-0.141, 'F':-0.023, 'P':0.286, 'S':0.238, 'T':0.147, 'W':-0.186, 'Y':-0.057, 'V':0.136}
    return length[aa_code]

def get_aa_alphapropensity(aa_code):
    alpha = {'A':-0.475, 'R':0.107, 'N':0.117, 'D':-0.014, 'C':0.070, 'Q':-0.030, 'E':-0.315, 'G':-0.024, 'H':0.041, 'I':-0.088, 'L':-0.265, 'K':-0.044, 'M':-0.155, 'F':0.072, 'P':0.407, 'S':-0.015, 'T':-0.015, 'W':0.389, 'Y':0.425, 'V':-0.187}
    return alpha[aa_code]

def get_aa_numbercodons(aa_code):
    ncodons = {'A':-0.039, 'R':-0.258, 'N':0.118, 'D':0.225, 'C':0.565, 'Q':0.035, 'E':0.157, 'G':0.018, 'H':0.280, 'I':-0.195, 'L':-0.274, 'K':-0.325, 'M':0.321, 'F':-0.002, 'P':-0.215, 'S':-0.068, 'T':-0.132, 'W':0.083, 'Y':-0.096, 'V':-0.196}
    return n_codons[aa_code]

def get_aa_betapropensity(aa_code):
    beta = {'A':0.181, 'R':-0.364, 'N':-0.055, 'D':0.156, 'C':-0.374, 'Q':-0.112, 'E':0.303, 'G':0.106, 'H':-0.021, 'I':-0.107, 'L':0.206, 'K':-0.027, 'M':0.077, 'F':0.208, 'P':0.384, 'S':-0.196, 'T':-0.274, 'W':0.297, 'Y':-0.091, 'V':-0.299}
    return beta[aa_code]

