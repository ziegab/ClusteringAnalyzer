import ROOT
from ROOT import TH2F, TCanvas
import numpy as np

def find_top_two_bins(hist, use_weighted=False):
    max_bin1 = (-1, -1)  # (x, y) coordinates of the highest bin
    max_bin2 = (-1, -1)  # (x, y) coordinates of the second highest bin
    max_content1 = -1
    max_content2 = -1
    
    nBinsX = hist.GetNbinsX()
    nBinsY = hist.GetNbinsY()

    for i in range(1, nBinsX + 1):
        for j in range(1, nBinsY + 1):
            if use_weighted:
                # Sum over neighboring bins (3x3 kernel)
                content = sum(hist.GetBinContent(i + dx, j + dy)
                              for dx in [-1, 0, 1] for dy in [-1, 0, 1]
                              if 1 <= i + dx <= nBinsX and 1 <= j + dy <= nBinsY)
            else:
                content = hist.GetBinContent(i, j)

            if content > max_content1:
                # Shift first max to second max
                max_content2 = max_content1
                max_bin2 = max_bin1

                # Update first max
                max_content1 = content
                max_bin1 = (i, j)
            elif content > max_content2:
                # Update second max
                max_content2 = content
                max_bin2 = (i, j)

    return [max_bin1, max_bin2]

def rotate_TH2(hist, bin_center=(8, 8), angle=None):
    """ Rotates a TH2 histogram around a given bin center.
        If angle is None, it is computed based on the two highest bins.
    """
    nBinsX = hist.GetNbinsX()
    nBinsY = hist.GetNbinsY()

    # Get the bin center coordinates
    x_center = hist.GetXaxis().GetBinCenter(bin_center[0])
    y_center = hist.GetYaxis().GetBinCenter(bin_center[1])

    # If angle is not provided, compute it from the top two bins
    if angle is None:
        top_bins = find_top_two_bins(hist)
        x1, y1 = hist.GetXaxis().GetBinCenter(top_bins[0][0]), hist.GetYaxis().GetBinCenter(top_bins[0][1])
        x2, y2 = hist.GetXaxis().GetBinCenter(top_bins[1][0]), hist.GetYaxis().GetBinCenter(top_bins[1][1])

        # Compute rotation angle to align (x1,y1) and (x2,y2) horizontally
        delta_y = y2 - y1
        delta_x = x2 - x1
        angle = np.arctan2(delta_y, delta_x)  # Angle in radians

    # print(f"Rotating by {np.degrees(angle):.2f} degrees")

    # Create new histogram with the same binning
    hRot = ROOT.TH2D("hRot", "Rotated Histogram", nBinsX, hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax(),
                     nBinsY, hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax())

    # Rotate each bin content
    for i in range(1, nBinsX + 1):
        for j in range(1, nBinsY + 1):
            x = hist.GetXaxis().GetBinCenter(i)
            y = hist.GetYaxis().GetBinCenter(j)
            content = hist.GetBinContent(i, j)

            # Translate to rotation center
            x_shifted = x - x_center
            y_shifted = y - y_center

            # Apply rotation
            x_rot = x_shifted * np.cos(-angle) - y_shifted * np.sin(-angle)
            y_rot = x_shifted * np.sin(-angle) + y_shifted * np.cos(-angle)

            # Translate back
            x_final = x_rot + x_center
            y_final = y_rot + y_center

            # Find nearest bin in rotated histogram
            binX_new = hRot.GetXaxis().FindBin(x_final)
            binY_new = hRot.GetYaxis().FindBin(y_final)

            if 1 <= binX_new <= nBinsX and 1 <= binY_new <= nBinsY:
                hRot.SetBinContent(binX_new, binY_new, content)

    return hRot

# # Get mass of each file from the csv file name
# def get_mass(csvfile, header):
#     posind = csvfile.find(header) + len(header) + 11
#     if posind != -1:
#         remaining_filename = csvfile[posind:].strip()
#         pmass = remaining_filename.split('m')[0]
#         mass = pmass.replace('p', '.')
#         return mass
#     else:
#         return None

def makeentry(l, output):
    index = int(l[0])
    gam = (float(l[1]))
    E = (float(l[2]))
    eta = float(l[3])
    phi = float(l[4])
    # deltaR0 = float(l[5])
    # deltaR1 = float(l[6])
    # C_2 = float(l[7])
    n_bins = 32
    length = float(int(n_bins/2))
    mH = E/gam
    H = TH2F("H", "#gamma = "+str(gam)+", E = "+str(E)+ " GeV, #eta = "+"{:.2f}".format(round(eta, 2))+", #phi = "+"{:.2f}".format(round(phi, 2))+";#eta;#phi", n_bins,-length,length,n_bins,-length,length)
    H.SetStats(0)
    for i in range(5, len(l), 3):
        H.Fill(int(l[i+1]),int(l[i+2]),float(l[i]))
    # maxcrystalE = H.GetBinContent(8, 8)
    # H.Scale(1.0 / maxcrystalE)
    H.Scale(1.0 / E)
    # H.Scale(1.0 / 1001.0)
    H_rotated = rotate_TH2(H, bin_center=(int(length)+1, int(length)+1), angle=None)
    array = H.GetArray()
    scaleeta = (eta+3)/6.0
    scalephi = (phi+np.pi)/(2*np.pi)
    with open(output+'.csv', mode='a') as file:
        file.write(str(index))
        file.write(", " + str(gam))
        file.write(", " + str(mH))
        file.write(", " + str(E))
        file.write(", " + str(scaleeta))
        file.write(", " + str(scalephi))
        # file.write(", " + str(C_2))
        for y in range(n_bins,0,-1):
            # for x in range(17,0,-1):
            for x in range(1, (n_bins+1)):
                index = x  + y * (n_bins+2)
                file.write(", " + str(array[index]))
                # print(f"{array[index]:.2f}", end=" ")
        file.write("\n")
        # print("")

import csv
import sys
import glob

file_dir = str(sys.argv[1])
print(file_dir)
csv_files = glob.glob(f"{file_dir}/*.csv")
file_counter = 0

for arg in sorted(csv_files):
    file_counter += 1
    startidx1 = arg.find("v")
    startidx2 = arg.find("v", startidx1 + 1)
    endidx = arg.rfind(".csv")
    versionname = arg[(startidx2 + 1):endidx]
    with open(arg) as csvfile:
        reader = csv.reader(csvfile)
        # mH = float(get_mass(str(sys.argv[1]), "SAIpreproc_AtoGG_"))
        output = "SAI_AtoGG_" + sys.argv[2] + "_MoE_32_Energy_totE_v" + str(versionname)
        n = 0
        for row in reader:
            n+=1
            # makeentry(row, output, mH)
            makeentry(row, output)

# with open(sys.argv[1]) as csvfile:
#         reader = csv.reader(csvfile)
#         # mH = float(get_mass(str(sys.argv[1]), "SAIpreproc_AtoGG_"))
#         output = sys.argv[2]
#         n = 0
#         for row in reader:
#             n+=1
#             # makeentry(row, output, mH)
#             makeentry(row, output)