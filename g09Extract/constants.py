import numpy as np
# 6-31G basis
basissize_631g = {'A':100, 'T':93, 'pPt':262, 'porphyrin':244, 'ethylene': 26, 'alkane4':87, 'alkane5':100, 'alkane6':113, 'alkane7':126, 'alkane8':139, 'alkane9':152, 'alkane10':165, 'alkane11':178,
             'alkane12':191, 'alkane13':204, "poly_repeats_2_orth": 145, "poly_repeats_3_orth": 202}

basissize_sto3g = {'A':55, "G":60, "T":51, "C":45, "3wj":1353, "porphyrin": 134}

basissize_zindo = {'A':45, 'T': 42, 'G': 49, 'C': 37, 'porphyrin': 110}
electronsize_zindo = {"A": 25, "T": 24, "G": 28, "C": 21, "porphyrin": 57}

electronsize_631g={'A':35, 'T':33, 'pPt':89, 'porphyrin':81, 'ethylene': 8}
electronsize_sto3g = {'A':35, 'G':39, 'T':33, 'C': 29, 'porphyrin':81}

# homos, lumos in order of D, B1, B2, ... A. zero-based
# b4_HOMOs = np.array([26, 10, 11, 12, 9, 25])
# b4_LUMOs = np.array([27, 31, 30, 29, 32, 28])

b4_symm_HOMOs = np.array([26, 9, 11, 12, 10, 25])
b4_symm_LUMOs = np.array([28, 32, 29, 30, 31, 27])

b5_symm_HOMOs = np.array([30, 11, 13, 12, 14, 10, 29])
b5_symm_LUMOs = np.array([32, 37, 34, 35, 33, 36, 31])

b6_symm_HOMOs = np.array([33, 12, 15, 14, 13, 16, 11, 34])
b6_symm_LUMOs = np.array([35, 41, 37, 40, 39, 38, 42, 36])

# b7_HOMOs = np.array([38, 12, 18, 14, 16, 15, 17, 13, 37])
# b7_LUMOs = np.array([40, 47, 41, 45, 43, 44, 42, 46, 39])

b7_symm_HOMOs=np.array([38, 13, 18, 14, 16, 15, 17, 12, 37])
b7_symm_LUMOs=np.array([40, 46, 42, 44, 43, 45, 41, 47, 39])

b7_HOMOs1 = np.array([37, 13, 18, 15, 16, 17, 14, 12, 38])
b7_LUMOs1 = np.array([40, 46, 41, 45, 43, 44, 42, 47, 39])

b7_HOMOs2 = np.array([37, 13, 14, 15, 16, 17, 18, 12, 38])
b7_LUMOs2 = np.array([39, 46, 45, 44, 43, 42, 41, 47, 40])

b7_HOMOs3 = np.array([38, 13, 18, 16, 15, 14, 17, 12, 37])
b7_LUMOs3 = np.array([39, 46, 41, 43, 44, 45, 42, 47, 40])

b7_HOMOs4 = np.array([37, 12, 17, 14, 15, 16, 18, 13, 38])
b7_LUMOs4 = np.array([39, 47, 41, 45, 43, 44, 42, 46, 40])

b7_HOMOs5 = np.array([38, 13, 18, 16, 17, 15, 14, 12, 37])
b7_LUMOs5 = np.array([39, 47, 41, 45, 43, 44, 42, 46, 40])

b8_symm_HOMOs = np.array([41, 14, 19, 15, 17, 18, 16, 20, 13, 42])
b8_symm_LUMOs = np.array([43, 51, 45, 47, 49, 50, 48, 46, 52, 44])

b9_symm_HOMOs = np.array([45, 15, 22, 17, 19, 20, 18, 16, 21, 14, 46])
b9_symm_LUMOs = np.array([48, 56, 50, 52, 54, 55, 53, 51, 49, 57, 47])

b10_symm_HOMOs = np.array([49, 16, 24, 17, 20, 21, 22, 19, 18, 23, 15, 50])
b10_symm_LUMOs = np.array([51, 61, 53, 55, 57, 59, 60, 58, 56, 54, 62, 52])

b11_symm_HOMOs = np.array([54, 16, 26, 21, 18, 23, 22, 24, 19, 20, 25, 17, 53])
b11_symm_LUMOs = np.array([55, 67, 57, 60, 64, 61, 65, 62, 63, 59, 58, 66, 56])

b12_symm_HOMOs = np.array([57, 18, 27, 19, 21, 23, 25, 26, 24, 22, 20, 28, 17, 58])
b12_symm_LUMOs = np.array([59, 72, 61, 63, 65, 67, 69, 70, 68, 66, 64, 62, 71, 60])

# b13_HOMOs = np.array([61, 19, 30, 20, 23, 24, 27, 28, 26, 25, 22, 21, 29, 18, 62])
# b13_LUMOs = np.array([64, 76, 65, 68, 69, 72, 73, 75, 74, 71, 70, 67, 66, 77, 63])

b13_symm_HOMOs=np.array([61, 18, 29, 20, 22, 24, 26, 28, 27, 25, 23, 21, 30, 19, 62])
b13_symm_LUMOs=np.array([64, 77, 65, 67, 69, 71, 73, 75, 74, 72, 70, 68, 66, 76, 63])

b13_HOMOs1 = np.array([62, 19, 24, 21, 26, 23, 29, 22, 25, 28, 30, 27, 20, 18, 61])
b13_LUMOs1 = np.array([64, 77, 65, 75, 67, 73, 68, 74, 71, 70, 69, 66, 72, 76, 63])

b13_HOMOs2 = np.array([62, 19, 24, 26, 23, 21, 27, 22, 28, 25, 29, 20, 30, 18, 61])
b13_LUMOs2 = np.array([63, 77, 67, 71, 72, 75, 69, 74, 70, 73, 66, 68, 65, 76, 64])

b13_HOMOs3 = np.array([61, 19, 29, 21, 27, 30, 28, 26, 24, 25, 22, 23, 20, 18, 62])
b13_LUMOs3 = np.array([64, 76, 65, 67, 70, 66, 68, 74, 72, 75, 71, 73, 69, 77, 63])

b13_HOMOs4 = np.array([61, 18, 24, 20, 23, 30, 28, 25, 22, 27, 26, 21, 29, 19, 62])
b13_LUMOs4 = np.array([63, 76, 65, 73, 69, 67, 70, 75, 71, 74, 68, 72, 66, 77, 64])

b13_HOMOs5 = np.array([61, 18, 27, 26, 28, 29, 30, 25, 22, 24, 21, 20, 23, 19, 62])
b13_LUMOs5 = np.array([64, 76, 65, 69, 72, 68, 67, 75, 71, 74, 70, 73, 66, 77, 63])

polynorbornyl_repeats_2_HOMOs=np.array([42, 14, 21, 19, 17, 15, 20, 18, 16, 22, 25, 24, 23, 26, 41])
polynorbornyl_repeats_2_LUMOs=np.array([43, 56, 50, 51, 55, 54, 52, 53, 57, 49, 47, 46, 48, 45, 44])

polynorbornyl_repeats_3_HOMOs = np.array([59, 22, 28, 26, 23, 27, 20, 21, 29, 25, 24, 30, 19, 34, 33, 31, 32, 36, 35, 37, 38, 60])
polynorbornyl_repeats_3_LUMOs = np.array([62, 81, 72, 76, 78, 70, 79, 82, 71, 75, 77, 69, 80, 68, 67, 74, 73, 65, 66, 64, 63, 61])

polynorbornyl_repeats_2_HOMOs_orth = np.array([45, 18, 22, 19, 25, 33, 17, 21, 20, 26, 27, 24, 23, 32, 31, 30, 44])
polynorbornyl_repeats_2_LUMOs_orth = np.array([47, 61, 57, 60, 48, 53, 62, 58, 59, 49, 56, 55, 54, 50, 51, 52, 46])


