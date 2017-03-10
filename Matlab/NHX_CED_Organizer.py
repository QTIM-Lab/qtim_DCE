import glob
import fnmatch
import os
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from shutil import copy, move
# import nifti_util
import csv

def OrganizeDSC():

    # outputpath = "C:/Users/azb22/Documents/Scripting/NHX/DSC_BATCH/"
    outputpath = "C:/Users/azb22/Documents/Scripting/CED/DSC_BATCH/"
    # inputpath = "C:/Users/azb22/Documents/Scripting/NHX/NORDIC_ICE"
    inputpath = "C:/Users/azb22/Documents/Scripting/CED/NORDIC_ICE"

    matches = []
    for root, dirnames, filenames in os.walk(inputpath):
        for filename in fnmatch.filter(filenames, 'ge_dsc_mc.nii'):
            matches.append(os.path.join(root, filename))

    for filename in matches:
        if "VISIT_01" in filename or "VISIT_02" in filename:
            print filename
            splitfile = str.split(filename, '\\')
            # new_filename = "_".join(splitfile[-3:])
            print inputpath + "GE_DSC_NII" + '/' + splitfile[-1]
            try:
                os.mkdir(inputpath + '\\' + "\\".join(splitfile[-3:-1]) + "\\GE_DSC_NII")
            except:
                pass
            copy(filename, inputpath + '\\' +  "\\".join(splitfile[-3:-1]) + "\\GE_DSC_NII" + '/' + splitfile[-1])

def CopyAutoAIFs():

    matches = []
    for root, dirnames, filenames in os.walk('/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/'):
        for filename in fnmatch.filter(filenames, 'dce_mc_st_eco1.nii'):
            matches.append(os.path.join(root, filename))

    for match in matches:
        if 'VISIT_01' in match or 'VISIT_02' in match:
            print match
            splitfile = str.split(match, '/')
            output_directory = '/'.join(splitfile[0:-1])
            visit = splitfile[-3]
            patient = splitfile[-4]
            aif_directory = '/qtim/users/data/CED/ANALYSIS/DCE/DCE_diff_AIFs/' + patient
            aiffile = visit + '_autoAIF_bAIF.txt'
            try:
                copy(aif_directory + '/' + aiffile, output_directory + '/' + aiffile)
            except:
                print 'ERROR'

    matches = []
    for root, dirnames, filenames in os.walk('/qtim2/users/data/NHX/ANALYSIS/DCE'):
        for filename in fnmatch.filter(filenames, 'dce_mc_st_eco1.nii'):
            matches.append(os.path.join(root, filename))

    for match in matches:
        if 'VISIT_01' in match or 'VISIT_02' in match:
            print match
            splitfile = str.split(match, '/')
            output_directory = '/'.join(splitfile[0:-1])
            visit = splitfile[-3]
            patient = splitfile[-4]
            aif_directory = '/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_diff_AIFs/NHX/' + patient
            aiffile = visit + '_autoAIF_bAIF.txt'
            try:
                copy(aif_directory + '/' + aiffile, output_directory + '/' + aiffile)
            except:
                print 'ERROR'

def CopyDCEs():
    destination_directory =  '/home/abeers/Projects/DCE_Method_Comparison/DCE_Aggregate'

    matches = []
    for root, dirnames, filenames in os.walk('/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/'):
        for filename in fnmatch.filter(filenames, 'dce_mc_st_eco1.nii'):
            matches.append(os.path.join(root, filename))

    for match in matches:
        if 'VISIT_01' in match or 'VISIT_02' in match:
            print match
            splitfile = str.split(match, '/')
            output_directory = '/'.join(splitfile[0:-1])
            visit = splitfile[-3]
            patient = splitfile[-4]
            try:
                copy(match, destination_directory + '/' + patient + '_' + visit + '_' + splitfile[-1])
            except:
                print 'ERROR'

    matches = []
    for root, dirnames, filenames in os.walk('/qtim2/users/data/NHX/ANALYSIS/DCE'):
        for filename in fnmatch.filter(filenames, 'dce_mc_st_eco1.nii'):
            matches.append(os.path.join(root, filename))

    for match in matches:
        if 'VISIT_01' in match or 'VISIT_02' in match:
            print match
            splitfile = str.split(match, '/')
            output_directory = '/'.join(splitfile[0:-1])
            visit = splitfile[-3]
            patient = splitfile[-4]
            try:
                copy(match, destination_directory + '/' + patient + '_' + visit + '_' + splitfile[-1])
            except:
                print 'ERROR'

def ROI_Scatter():

    path = '/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/ROI_ONLY/'

    # ROI_list = glob.glob(path + '*')
    ROI_list = []
    for root, dirnames, filenames in os.walk(path):
        for filename in filenames:
            ROI_list.append(os.path.join(root, filename))

    for ROI_idx, ROI in enumerate(ROI_list):
        if 'ktrans' in ROI:

            print ROI

            ktrans_nifti = nib.load(ROI)
            ktrans_numpy = ktrans_nifti.get_data()

            filename = str.split(ROI, '/')[-1]
            ve_ROI = ROI[0:-13] + 've.nii.gz'

            ve_nifti = nib.load(ve_ROI)
            ve_numpy = ve_nifti.get_data()

            ve_flatten = ve_numpy.flatten()
            ktrans_flatten = ktrans_numpy.flatten()

            plt.scatter(ktrans_flatten, ve_flatten)
            plt.xlim(0.01, .5)
            plt.ylim(0.01, .5)
            plt.show()

def ROI_Repeatability_Package(regex, param):

    path = '/home/abeers/Projects/DCE_Package/Matlab/Test_Results/Corrected/'

    # regex = 'High_Blur'

    ROI_list = sorted(glob.glob(path + '*' + regex + '*' + param + '*'))
    # ROI_list = []
    # for root, dirnames, filenames in os.walk(path):
        # for filename in filenames:
            # ROI_list.append(os.path.join(root, filename))

    visits = ['visit_01', 'visit_02']
    values = ['mean', 'median', 'min', 'max','std']
    
    statistics_worksheet = np.zeros((len(ROI_list) + 1,len(values)+1), dtype=object)  
    statistics_worksheet[0,:] = ['name'] + values

    for ROI_idx, ROI in enumerate(ROI_list):
        # print ROI
        ROI_nifti = nib.load(ROI)
        ROI_numpy = ROI_nifti.get_data()

        filename = str.split(ROI, '/')[-1]
        patient = filename[0:6]
        visit = filename[7:15]
        ID = patient + '_' + visit

        masked_ROI_numpy = np.ma.masked_where(ROI_numpy <= 0, ROI_numpy)
        masked_ROI_numpy = np.ma.masked_where(masked_ROI_numpy >= 1, masked_ROI_numpy)

        ROI_values = [np.ma.mean(masked_ROI_numpy), np.ma.median(masked_ROI_numpy), np.ma.min(masked_ROI_numpy), np.ma.max(masked_ROI_numpy), np.ma.std(masked_ROI_numpy)]
        # print ROI_values

        statistics_worksheet[ROI_idx+1,0] = ID
        statistics_worksheet[ROI_idx+1, 1:6] = ROI_values

    print statistics_worksheet
    row_num = 0
    visit_1_means = []
    visit_2_means = []
    visit_1_medians = []
    visit_2_medians = []
    for row_idx, row in enumerate(statistics_worksheet):
        if row_idx != 0:
            print row
            if 'VISIT_01' in row[0] and not np.ma.is_masked(row[1]):
                if 'VISIT_02' in statistics_worksheet[row_idx+1, 0] and row[0][0:7] in statistics_worksheet[row_idx+1,0]:
                    visit_1_means += [float(row[1])]
                    visit_2_means += [float(statistics_worksheet[row_idx+1,1])]
                    visit_1_medians += [float(row[2])]
                    visit_2_medians += [float(statistics_worksheet[row_idx+1,2])]

    mean_CR = CalcCR(visit_1_means, visit_2_means)
    median_CR = CalcCR(visit_1_medians, visit_2_medians)

    with open(param + regex + '_ROI_statistics_' + str(mean_CR) + '_' + str(median_CR) + '.csv', 'wb') as writefile:
        csvfile = csv.writer(writefile, delimiter=',')
        for row in statistics_worksheet:
            if str(row[0]) != '0':
                csvfile.writerow(row)












def ROIStatistics(param='ktrans'):

    path = '/qtim2/users/data/NHX/ANALYSIS/DCE/DCE_KTRANS_ANALYSIS/PCA_BLUR_TESTING/'

    # ROI_list = glob.glob(path + '*')
    ROI_list = []
    for root, dirnames, filenames in os.walk(path):
        for filename in filenames:
            ROI_list.append(os.path.join(root, filename))

    echos = ['Corrected']
    visits = ['visit_01', 'visit_02']
    kernels = ['0','.65','1.35']
    PCAs = ['0','5','10']
    headers = ['static_populationAIF','T1Map_populationAIF','static_automaticAIF','T1MAP_automaticAIF']
    values = ['mean', 'median', 'min', 'max','std']
    
    statistics_worksheet = np.zeros((2000,1+len(echos)*len(kernels)*len(PCAs)*len(headers)*len(values)*len(visits)), dtype=object)

    statistics_worksheet[0,0] = 'Patient_ID'
    for echo_idx, echo in enumerate(echos):
        for visit_idx, visit in enumerate(visits):
            for kernel_idx, kernel in enumerate(kernels):
                for PCA_idx, PCA in enumerate(PCAs):
                    for header_idx, header in enumerate(headers):
                        for value_idx, value in enumerate(values):
                            statistics_worksheet[0, 1 + echo_idx*len(visits)*len(kernels)*len(PCAs)*len(headers)*len(values) + visit_idx*len(kernels)*len(PCAs)*len(headers)*len(values)+kernel_idx*len(PCAs)*len(headers)*len(values) + PCA_idx*len(values)*len(headers)+len(values)*header_idx + value_idx] = echo + '_' + visit + '_' + kernel + '_' + PCA + '_' + header + '_' + value

    for ROI_idx, ROI in enumerate(ROI_list):
        if param in ROI:
            print ROI
            ROI_nifti = nib.load(ROI)
            ROI_numpy = ROI_nifti.get_data()

            filename = str.split(ROI, '/')[-1]
            patient = filename[0:6]
            visit = filename[7:15]
            ID = patient + '_' + visit

            masked_ROI_numpy = np.ma.masked_where(ROI_numpy <= 0, ROI_numpy)
            # masked_ROI_numpy = np.ma.masked_where(masked_ROI_numpy == 3, masked_ROI_numpy)

            ROI_values = [np.ma.mean(masked_ROI_numpy), np.ma.median(masked_ROI_numpy), np.ma.min(masked_ROI_numpy), np.ma.max(masked_ROI_numpy), np.ma.std(masked_ROI_numpy)]
            print ROI_values

            colid = 1


            if 'autoAIF' in ROI:
                colid += 10

            if 't1map' in ROI:
                colid += 5

            if 'PCA_5' in ROI:
                colid += 20
            elif 'PCA_10' in ROI:
                colid += 40

            if 'kernel_0.65' in ROI:
                colid += 60
            if 'kernel_1.35' in ROI:
                colid += 120

            if 'VISIT_02' in ROI:
                colid += 180

            # if 'Echo1' in ROI:
            #     colid += 40
            # elif 'Echo2' in ROI:
            #     colid += 80

            current_rowid = 1
            for row_idx, row in enumerate(statistics_worksheet):
                if str(row[0]) == '0' or str(row[0]) == patient:
                    current_rowid = row_idx
                    break

            statistics_worksheet[current_rowid,0] = patient
            statistics_worksheet[current_rowid, colid:colid+5] = ROI_values

    emptied_statistics_worksheet = np.copy(statistics_worksheet)
    emptied_statistics_worksheet[statistics_worksheet==0] == ''

    with open(param + '_method_statistics.csv', 'wb') as writefile:
        csvfile = csv.writer(writefile, delimiter=',')
        for row in emptied_statistics_worksheet:
            if str(row[0]) != '0':
                csvfile.writerow(row)

    new_statistics_worksheet = np.copy(statistics_worksheet)
    test_masked_statistics_worksheet = np.ma.masked_where(new_statistics_worksheet==0, new_statistics_worksheet)
    test_masked_statistics_worksheet = np.ma.masked_where(new_statistics_worksheet>1, test_masked_statistics_worksheet)
    with open(param + '_mask_test.csv', 'wb') as writefile:
        csvfile = csv.writer(writefile, delimiter=',')
        for row in test_masked_statistics_worksheet:
            if str(row[0]) != '0':
                csvfile.writerow(row)

    masked_statistics_worksheet = test_masked_statistics_worksheet

    combined_statistics_worksheet = np.copy(statistics_worksheet)

    # for row_idx, row in enumerate(statistics_worksheet):
    #     if str(row[0]) == '0':
    #         combined_statistics_worksheet[row_idx:row_idx*2-1,1:1+len(kernels)*len(PCAs)*len(headers)*len(values)*len(visits)/2] = statistics_worksheet[1:row_idx,1+len(kernels)*len(PCAs)*len(headers)*len(values)*len(visits)/2:]
    #         break

    # CCC_columns = []
    # for col in xrange(statistics_worksheet.shape[1]):
    #     # if (col-1) % 5 < 2 and col > 0:
    #     if (col-1) % len(values) == 1 and col > 0 and col < int(statistics_worksheet.shape[1]/2):
    #         CCC_columns += [col]
    # CCC_statistics = np.zeros((len(CCC_columns) + 1, len(CCC_columns) + 1), dtype=object)

    # CCC_statistics[0,0] = 'CCC_Statistic'

    # for stat1_idx, stat1 in enumerate(CCC_columns):
    #     print statistics_worksheet[0,stat1]
    #     CCC_statistics[stat1_idx+1,0] = statistics_worksheet[0,stat1]
    #     CCC_statistics[0,stat1_idx+1] = statistics_worksheet[0,stat1]

    # for stat1_idx, stat1 in enumerate(CCC_columns):
    #     for stat2_idx, stat2 in enumerate(CCC_columns):
    #         CCC_statistics[stat1_idx + 1, stat2_idx + 1] = CalcCCC(masked_statistics_worksheet[1:,stat1], masked_statistics_worksheet[1:,stat2])

    CR_columns = []
    for col in xrange(statistics_worksheet.shape[1]):
        if (col-1) % len(values) < 2 and col > 0:
            CR_columns += [col]
    CR_statistics = np.zeros((len(CR_columns) + 1, len(CR_columns) + 1), dtype=object)

    CR_statistics[0,0] = 'CR_Statistic'

    for stat1_idx, stat1 in enumerate(CR_columns):
        print statistics_worksheet[0,stat1]
        CR_statistics[stat1_idx+1,0] = statistics_worksheet[0,stat1]
        CR_statistics[0,stat1_idx+1] = statistics_worksheet[0,stat1]

    for stat1_idx, stat1 in enumerate(CR_columns):
        for stat2_idx, stat2 in enumerate(CR_columns):
            CR_statistics[stat1_idx + 1, stat2_idx + 1] = CalcCR(masked_statistics_worksheet[1:,stat1], masked_statistics_worksheet[1:,stat2])
            print [stat1, stat2]
            print masked_statistics_worksheet[1:,stat1]
            print masked_statistics_worksheet[1:,stat2]
            print 'next..'

    # with open(param + '_CCC_statistics.csv', 'wb') as writefile:
    #     csvfile = csv.writer(writefile, delimiter=',')
    #     for row in CCC_statistics:
    #         csvfile.writerow(row)

    with open(param + '_CR_statistics.csv', 'wb') as writefile:
        csvfile = csv.writer(writefile, delimiter=',')
        for row in CR_statistics:
            csvfile.writerow(row)

    return

def T1_CR(param="t1"):

    matches = []

    for vol_roots in ['/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/', '/qtim2/users/data/NHX/ANALYSIS/DCE/']:
        for root, dirnames, filenames in os.walk(vol_roots):
            for filename in fnmatch.filter(filenames, 'rT1AxialPostROI.nii'):
                if 'VISIT_01' in os.path.join(root, filename) or 'VISIT_02' in os.path.join(root, filename):
                    matches.append(os.path.join(root, filename))

    ROI_list = matches

    visits = ['visit_01', 'visit_02']
    values = ['mean', 'median', 'min', 'max','std']
    
    statistics_worksheet = np.zeros((2000,1+len(values)*len(visits)), dtype=object)

    statistics_worksheet[0,0] = 'Patient_ID'
    for visit_idx, visit in enumerate(visits):
        for value_idx, value in enumerate(values):
            statistics_worksheet[0, 1+visit_idx*len(values) + value_idx] = visit + '_' + value

    for ROI_idx, ROI in enumerate(ROI_list):
        if 'rT1AxialPostROI.nii' in ROI:
            print ROI
            ROI_nifti = nib.load(ROI)
            ROI_numpy = ROI_nifti.get_data()

            splitfile = str.split(ROI, '/')

            filename = splitfile[-1]
            patient = splitfile[-5]
            visit = splitfile[-4]
            ID = patient + '_' + visit

            T1map = '/'.join(splitfile[0:-3]) + '/MAPS/T1inDCE.nii'

            try:
                t1_nifti = nib.load(T1map)
                t1_numpy = t1_nifti.get_data()
                t1_numpy = np.squeeze(t1_numpy)

                # print t1_numpy.shape
                # print ROI_numpy.shape
                # print np.unique(ROI_numpy)
                # print np.unique(t1_numpy)

                # nifti_util.check_image(t1_numpy, ROI_numpy, mode="maximal_slice")

                masked_ROI_numpy = np.ma.masked_where(ROI_numpy != 1, t1_numpy)

                ROI_values = [np.ma.mean(masked_ROI_numpy), np.ma.median(masked_ROI_numpy), np.ma.min(masked_ROI_numpy), np.ma.max(masked_ROI_numpy), np.ma.std(masked_ROI_numpy)]
                print ROI_values

                colid = 1

                if 'VISIT_02' in ROI:
                    colid += 5

                current_rowid = 1
                for row_idx, row in enumerate(statistics_worksheet):
                    if str(row[0]) == '0' or str(row[0]) == patient:
                        current_rowid = row_idx
                        break

                statistics_worksheet[current_rowid,0] = patient
                statistics_worksheet[current_rowid, colid:colid+5] = ROI_values
            except:
                pass

    emptied_statistics_worksheet = np.copy(statistics_worksheet)
    emptied_statistics_worksheet[statistics_worksheet==0] == ''

    with open(param + '_method_statistics.csv', 'wb') as writefile:
        csvfile = csv.writer(writefile, delimiter=',')
        for row in emptied_statistics_worksheet:
            csvfile.writerow(row)

    # combined_statistics_worksheet = np.copy(statistics_worksheet)

    # for row_idx, row in enumerate(statistics_worksheet):
    #     if str(row[0]) == '0':
    #         combined_statistics_worksheet[row_idx:row_idx*2-1,1:1+len(kernels)*len(PCAs)*len(headers)*len(values)*len(visits)/2] = statistics_worksheet[1:row_idx,1+len(kernels)*len(PCAs)*len(headers)*len(values)*len(visits)/2:]
    #         break

    # masked_statistics_worksheet = np.ma.masked_where(combined_statistics_worksheet>1, statistics_worksheet)
    # masked_statistics_worksheet = np.ma.masked_where(masked_statistics_worksheet==0, statistics_worksheet)

    # CCC_columns = []
    # for col in xrange(statistics_worksheet.shape[1]):
    #     # if (col-1) % 5 < 2 and col > 0:
    #     if (col-1) % 5 == 1 and col > 0 and col < int(statistics_worksheet.shape[1]/2):
    #         CCC_columns += [col]
    # CCC_statistics = np.zeros((len(CCC_columns) + 1, len(CCC_columns) + 1), dtype=object)

    # CCC_statistics[0,0] = 'CCC_Statistic'

    # for stat1_idx, stat1 in enumerate(CCC_columns):
    #     print statistics_worksheet[0,stat1]
    #     CCC_statistics[stat1_idx+1,0] = statistics_worksheet[0,stat1]
    #     CCC_statistics[0,stat1_idx+1] = statistics_worksheet[0,stat1]

    # for stat1_idx, stat1 in enumerate(CCC_columns):
    #     for stat2_idx, stat2 in enumerate(CCC_columns):
    #         CCC_statistics[stat1_idx + 1, stat2_idx + 1] = CalcCCC(masked_statistics_worksheet[1:,stat1], masked_statistics_worksheet[1:,stat2])

    statistics_worksheet = np.ma.masked_where(statistics_worksheet==0, statistics_worksheet)
    masked_statistics_worksheet = np.ma.masked_where(statistics_worksheet > 5000, statistics_worksheet)

    CR_columns = []
    for col in xrange(statistics_worksheet.shape[1]):
        if (col-1) % 5 == 1 and col > 0:
            CR_columns += [col]
    CR_statistics = np.zeros((len(CR_columns) + 1, len(CR_columns) + 1), dtype=object)

    CR_statistics[0,0] = 'CR_Statistic'

    for stat1_idx, stat1 in enumerate(CR_columns):
        print statistics_worksheet[0,stat1]
        CR_statistics[stat1_idx+1,0] = statistics_worksheet[0,stat1]
        CR_statistics[0,stat1_idx+1] = statistics_worksheet[0,stat1]

    for stat1_idx, stat1 in enumerate(CR_columns):
        for stat2_idx, stat2 in enumerate(CR_columns):
            CR_statistics[stat1_idx + 1, stat2_idx + 1] = CalcCR(masked_statistics_worksheet[1:,stat1], masked_statistics_worksheet[1:,stat2])

    # with open(param + '_CCC_statistics.csv', 'wb') as writefile:
    #     csvfile = csv.writer(writefile, delimiter=',')
    #     for row in CCC_statistics:
    #         csvfile.writerow(row)

    with open(param + '_CR_statistics.csv', 'wb') as writefile:
        csvfile = csv.writer(writefile, delimiter=',')
        for row in CR_statistics:
            csvfile.writerow(row)

    return

def MoveAutoAIFs():
    AIFS = np.genfromtxt('Auto_AIF_Comparison.csv',dtype='object',delimiter=',', missing_values='')
    print AIFS
    for col in xrange(AIFS.shape[1]):
        if col == 0:
            continue
        if AIFS[2,col] == '':
            continue
        code = AIFS[0,col]
        patient = code[0:2]
        visit = code[2:]
        AIF = AIFS[2:, col]
        output_filename = '/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/CED_' + str(patient) + '/VISIT_' + str(visit) + '/MAPS/NORDIC_ICE_AIF.txt'
        output_string = ''
        for item in AIF:
            output_string +=  item + ';'
        output_string = output_string[0:-1]
        print output_string
        f = open(output_filename,'wb')
        f.write(output_string)
        f.close()

def CalcCCC(x, y):
    mean_x = np.ma.mean(x)
    mean_y = np.ma.mean(y)
    std_x = np.ma.std(x)
    std_y = np.ma.std(y)
    correl = np.ma.corrcoef(x,y)[0,1]
    CCC = (2 * correl * std_x * std_y) / (np.ma.var(x) + np.ma.var(y) + np.square(mean_x - mean_y))
    return CCC

def CalcCR(x,y):
    z = np.zeros_like(x)
    for index in xrange(len(x)):
        z[index] = np.ma.power(x[index] - y[index], 2)

    # print z
    z = np.ma.sum(z) / len(z)
    z = np.ma.power(z, 0.5)
    z = 1.96 * z

    return z

def CollectK2Statistics():

    k2_statistics = np.zeros((2000,11),dtype=object)

    visits = ['visit_01', 'visit_02']
    values = ['mean', 'median', 'min', 'max','std']
    
    k2_statistics[0,0] = 'Patient_ID'
    for visit_idx, visit in enumerate(visits):
        for value_idx, value in enumerate(values):
            k2_statistics[0, 1+5*visit_idx+value_idx] = 'k2_' + visit + '_' + value


    matches = []
    for root, dirnames, filenames in os.walk('C:/Users/azb22/Documents/Scripting/CED/NORDIC_ICE'):
        for filename in fnmatch.filter(filenames, 'leakage_correct_params_Leakage*'):
            matches.append(os.path.join(root, filename))

    for match in matches:
        if 'VISIT_01' or 'VISIT_02' in match:
            print match
            k2_nifti = nib.load(match)
            k2_numpy = k2_nifti.get_data()

            k2_directory = '/'.join(str.split(match, '\\')[0:-1])

            try:
                ROI_nifti = nib.load(k2_directory + '/t1axialpostroi_lowres.nii')
                ROI_numpy = ROI_nifti.get_data()

                masked_ROI_numpy = np.ma.masked_where(ROI_numpy == 0, k2_numpy)
                masked_ROI_numpy = np.ma.masked_where(masked_ROI_numpy == 0, k2_numpy)

                ROI_values = [np.ma.mean(masked_ROI_numpy), np.ma.median(masked_ROI_numpy), np.ma.min(masked_ROI_numpy), np.ma.max(masked_ROI_numpy), np.ma.std(masked_ROI_numpy)]

                print ROI_values

                patient = 'CED_' + str.split(match, 'CED_')[-1][0:2]

                current_rowid = 1
                for row_idx, row in enumerate(k2_statistics):
                    if str(row[0]) == '0' or str(row[0]) == patient:
                        current_rowid = row_idx
                        break

                if 'VISIT_01' in match:
                    k2_statistics[current_rowid, 1:6] = ROI_values
                else:
                    k2_statistics[current_rowid, 6:] = ROI_values
            except:
                print 'ERROR! ' + patient

            current_rowid = 1
            for row_idx, row in enumerate(k2_statistics):
                if str(row[0]) == '0' or str(row[0]) == patient:
                    current_rowid = row_idx
                    break

            patient = 'CED_' + str.split(match, 'CED_')[-1][0:2]
            k2_statistics[current_rowid,0] = patient

    matches = []
    for root, dirnames, filenames in os.walk('C:/Users/azb22/Documents/Scripting/NHX/NORDIC_ICE'):
        for filename in fnmatch.filter(filenames, 'leakage_correct_params_Leakage*'):
            matches.append(os.path.join(root, filename))

    for match in matches:
        match_checker = '/'.join(str.split(match, '\\')[0:-1])
        if 'VISIT_01' or 'VISIT_02' in match_checker:
            print match
            k2_nifti = nib.load(match)
            k2_numpy = k2_nifti.get_data()

            k2_directory = '/'.join(str.split(match, '\\')[0:-1])

            try:
                ROI_nifti = nib.load(k2_directory + '/LowResT1AxialPostROI.nii')
                ROI_numpy = ROI_nifti.get_data()
                ROI_numpy = ROI_numpy[:,:,:,0]

                masked_ROI_numpy = np.ma.masked_where(ROI_numpy == 0, k2_numpy)
                masked_ROI_numpy = np.ma.masked_where(masked_ROI_numpy == 0, k2_numpy)

                ROI_values = [np.ma.mean(masked_ROI_numpy), np.ma.median(masked_ROI_numpy), np.ma.min(masked_ROI_numpy), np.ma.max(masked_ROI_numpy), np.ma.std(masked_ROI_numpy)]

                print ROI_values

                patient = 'NHX_' + str.split(match, 'NHX_')[-1][0:2]

                current_rowid = 1
                for row_idx, row in enumerate(k2_statistics):
                    if str(row[0]) == '0' or str(row[0]) == patient:
                        current_rowid = row_idx
                        break

                if 'VISIT_01' in match_checker:
                    k2_statistics[current_rowid, 1:6] = ROI_values
                else:
                    k2_statistics[current_rowid, 6:] = ROI_values
            except:
                print 'ERROR! ' + patient

            current_rowid = 1
            for row_idx, row in enumerate(k2_statistics):
                if str(row[0]) == '0' or str(row[0]) == patient:
                    current_rowid = row_idx
                    break

            patient = 'NHX_' + str.split(match, 'NHX_')[-1][0:2]
            k2_statistics[current_rowid,0] = patient

    k2_statistics[k2_statistics=='0'] == ''
    k2_statistics[k2_statistics=='--'] == ''


    with open('k2_statistics.csv', 'wb') as writefile:
        csvfile = csv.writer(writefile, delimiter=',')
        for row in k2_statistics:
            csvfile.writerow(row)

def r2star_creator(TE1=2.73, TE2=3.89):

    study_files = []

    for root, dirnames, filenames in os.walk('/qtim2/users/data/NHX/ANALYSIS/DCE/'):
        for filename in fnmatch.filter(filenames, 'dce_mc_st_eco1.nii'):
            study_files.append(os.path.join(root, filename))

    for root, dirnames, filenames in os.walk('/qtim/users/data/CED/ANALYSIS/DCE/PREPARATION_FILES/'):
        for filename in fnmatch.filter(filenames, 'dce_mc_st_eco1.nii'):
            study_files.append(os.path.join(root, filename))


    for study_file in study_files:
        if 'dce_mc_st_eco1.nii' in study_file and ('VISIT_01' in study_file or 'VISIT_02' in study_file):
            fileid = str.split(study_file, '/')

            study_file_2 = '/'.join(fileid[0:-1]) + '/dce_mc_st_eco2.nii'
            visit_id = fileid[-3]
            patient_id = fileid[-4]

            print [visit_id, patient_id]

            study_1_nifti = nib.load(study_file)
            study_2_nifti = nib.load(study_file_2)
            study_1_numpy = study_1_nifti.get_data()
            study_2_numpy = study_2_nifti.get_data()

            study_ratio = study_1_numpy / study_2_numpy
            study_ratio[np.isinf(study_ratio)] = 0

            R2star = np.log(study_ratio) / (TE2 - TE1)

            study_corrected = study_1_numpy * np.exp(TE1 * R2star)

            study_corrected_nifti = nib.Nifti1Image(study_corrected, study_1_nifti.get_affine())
            nib.save(study_corrected_nifti, '/'.join(fileid[0:-1]) + '/dce_mc_st_corrected.nii')





    return


if __name__ == "__main__":
    # r2star_creator()
    # CollectK2Statistics()
    # T1_CR()
    # ROI_Scatter()
    # ROIStatistics('ktrans')
    for regex in ['_Bad_Corrected']:
        for param in ['ktrans','auc','ve']:
            ROI_Repeatability_Package(regex, param)
    # ROI_Repeatability_Package('_No_Blur_','ktrans')
    # ROI_Repeatability_Package('_No_Blur_','ve')
    # ROI_Repeatability_Package('_No_Blur_','auc')
    # ROI_Repeatability_Package('_High_Blur_','ktrans')
    # ROI_Repeatability_Package('_High_Blur_','ve')
    # ROI_Repeatability_Package('_High_Blur_','auc')
    # ROI_Repeatability_Package('_Alt_Integrate_','ktrans')
    # ROI_Repeatability_Package('_Alt_Integrate_','ve')
    # ROI_Repeatability_Package('_Alt_Integrate_','auc')
    # ROIStatistics('ve')
    # MoveAutoAIFs()
    # CopyDCEs()
    # CopyAutoAIFs()
    # OrganizeDSC()