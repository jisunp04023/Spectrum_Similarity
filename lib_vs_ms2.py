import pandas as pd
import numpy as np
import math
from scipy import stats

lib = open("Lib_sample.txt", 'r')
lib_df = pd.DataFrame(columns=['seq', 'charge'])

while True:
    libline = lib.readline()

    if not libline: break

    if libline.startswith('SEQ='):
        seq = libline[4:-1]
        libline = lib.readline()
        charge = libline[7:-1]
        data = {'seq': [seq], 'charge': [charge]}
        data_df = pd.DataFrame(data)
        lib_df = lib_df.append(data_df)
lib.close()

lib_df = lib_df.reset_index(drop=True)

ms2 = open("MS2_sample.txt", 'r')
ms2_df = pd.DataFrame(columns=['seq', 'charge'])

while True:
    ms2line = ms2.readline()

    if not ms2line: break

    if ms2line.startswith('SEQ='):
        seq = ms2line[4:-1]
        ms2line = ms2.readline()
        ms2line = ms2.readline()
        charge = ms2line[7:-1]
        data = {'seq': [seq], 'charge': [charge]}
        data_df = pd.DataFrame(data)
        ms2_df = ms2_df.append(data_df)
ms2.close()

ms2_df = ms2_df.reset_index(drop=True)

isIn_df = pd.DataFrame(columns=['seq', 'charge', 'ms2_idx'])

for i in range(len(ms2_df)):
    # print(i)
    s = ms2_df.loc[i].values[0]
    ch = ms2_df.loc[i].values[1]

    isIn = lib_df[(lib_df['seq'] == s) & (lib_df['charge'] == ch)].copy()

    if len(isIn):
        # print(isIn)
        isIn['ms2_idx'] = i
        isIn_df = isIn_df.append(isIn)


def getMS2Spectrum(idx):
    ms2 = open("MS2_sample.txt", 'r')
    ms2spec_df = pd.DataFrame(columns=['m/z', 'intensity'])
    i = 0

    while True:
        ms2line = ms2.readline()

        if not ms2line: break

        if ms2line.startswith('SCANS='):
            if i == idx:

                while True:
                    ms2line = ms2.readline()
                    if ms2line.startswith('END'): break

                    data = ms2line[:-1].split('\t')
                    mz = float(ms2line[:-1].split('\t')[0])
                    inten = float(ms2line[:-1].split('\t')[1])
                    data = {'m/z': [mz], 'intensity': [math.sqrt(inten)]}
                    data_df = pd.DataFrame(data)
                    ms2spec_df = ms2spec_df.append(data_df)

            i += 1

    ms2.close()
    ms2spec_df = ms2spec_df.reset_index(drop=True)
    return ms2spec_df


def getLibSpectrum(idx):
    lib = open("Lib_sample.txt", 'r')
    libspec_df = pd.DataFrame(columns=['m/z', 'intensity'])
    i = 0

    while True:
        libline = lib.readline()

        if not libline: break

        if libline.startswith('CHARGE='):
            if i == idx:

                while True:
                    libline = lib.readline()
                    if libline.startswith('END'): break

                    data = libline[:-1].split('\t')
                    mz = float(libline[:-1].split('\t')[0])
                    inten = float(libline[:-1].split('\t')[1])
                    data = {'m/z': [mz], 'intensity': [math.sqrt(inten)]}
                    data_df = pd.DataFrame(data)
                    libspec_df = libspec_df.append(data_df)

            i += 1

    lib.close()
    libspec_df = libspec_df.reset_index(drop=True)
    return libspec_df


def PCC(MS2, LIB):
    mz = [0]
    lib = [0]
    ms2 = [0]

    libidx = 0
    ms2idx = 0

    while True:
        if libidx == len(LIB):
            while True:
                if ms2idx == len(MS2): break
                mz += [MS2[ms2idx][0]]
                lib += [0]
                ms2 += [MS2[ms2idx][1]]
                ms2idx += 1
            break
        if ms2idx == len(MS2):
            while True:
                if libidx == len(LIB): break
                mz += [LIB[libidx][0]]
                lib += [LIB[libidx][1]]
                ms2 += [0]
                libidx += 1
            break

        if abs(LIB[libidx][0] - MS2[ms2idx][0]) <= 0.05:
            mz += [(LIB[libidx][0] + MS2[ms2idx][0]) / 2]
            lib += [LIB[libidx][1]]
            ms2 += [MS2[ms2idx][1]]
            libidx += 1
            ms2idx += 1
        elif LIB[libidx][0] < MS2[ms2idx][0]:
            mz += [LIB[libidx][0]]
            lib += [LIB[libidx][1]]
            ms2 += [0]
            libidx += 1
        else:
            mz += [MS2[ms2idx][0]]
            lib += [0]
            ms2 += [MS2[ms2idx][1]]
            ms2idx += 1

    del mz[0]
    del lib[0]
    del ms2[0]

    return stats.pearsonr(lib, ms2)[0]


def DotProduct(MS2, LIB):
    mz = [0]
    lib = [0]
    ms2 = [0]

    libidx = 0
    ms2idx = 0

    while True:
        if libidx == len(LIB):
            while True:
                if ms2idx == len(MS2): break
                mz += [MS2[ms2idx][0]]
                lib += [0]
                ms2 += [MS2[ms2idx][1]]
                ms2idx += 1
            break
        if ms2idx == len(MS2):
            while True:
                if libidx == len(LIB): break
                mz += [LIB[libidx][0]]
                lib += [LIB[libidx][1]]
                ms2 += [0]
                libidx += 1
            break

        if abs(LIB[libidx][0] - MS2[ms2idx][0]) <= 0.05:
            mz += [(LIB[libidx][0] + MS2[ms2idx][0]) / 2]
            lib += [LIB[libidx][1]]
            ms2 += [MS2[ms2idx][1]]
            libidx += 1
            ms2idx += 1
        elif LIB[libidx][0] < MS2[ms2idx][0]:
            mz += [LIB[libidx][0]]
            lib += [LIB[libidx][1]]
            ms2 += [0]
            libidx += 1
        else:
            mz += [MS2[ms2idx][0]]
            lib += [0]
            ms2 += [MS2[ms2idx][1]]
            ms2idx += 1

    del mz[0]
    del lib[0]
    del ms2[0]

    similarity = np.dot(lib, ms2) / math.sqrt(np.dot(lib, lib) * np.dot(ms2, ms2))

    return similarity


similarity_df = pd.DataFrame(columns=['seq', 'charge', 'PCC', 'DotProduct'])

for i in isIn_df.index:
    s = isIn_df.loc[i, 'seq']
    ch = isIn_df.loc[i, 'charge']

    j = isIn_df.loc[i, 'ms2_idx']

    MS2 = getMS2Spectrum(j).values.tolist()
    LIB = getLibSpectrum(i).values.tolist()

    pcc = PCC(MS2, LIB)
    dot = DotProduct(MS2, LIB)

    data = {'seq': [s], 'charge': [ch], 'PCC': [pcc], 'DotProduct': [dot]}
    data_df = pd.DataFrame(data)
    similarity_df = similarity_df.append(data_df)

similarity_df.to_excel("Similarity.xlsx")