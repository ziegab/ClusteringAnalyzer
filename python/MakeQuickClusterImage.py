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

def makeplot(l, n):
    gam = int(float(l[0]))
    E = int(float(l[1]))
    eta = float(l[2])
    phi = float(l[3])
    # mass = float(l[1])/float(l[0])
    n_bins = 32
    length = float(int(n_bins/2))
    H = TH2F("H", "#gamma = "+str(gam)+", E = "+str(E)+ " GeV, #eta = "+"{:.2f}".format(round(eta, 2))+", #phi = "+"{:.2f}".format(round(phi, 2))+";#eta;#phi", n_bins,-length,length,n_bins,-length,length)
    H.SetStats(0)
    for i in range(4, len(l), 3):
        H.Fill(int(l[i+1]),int(l[i+2]),float(l[i]))
    # maxcrystalE = H.GetBinContent(8, 8)
    if float(l[1]) > 0:
        H.Scale(1.0 / float(l[1]))

    H_rotated = rotate_TH2(H, bin_center=(int(length)+1, int(length)+1), angle=None)

    C = TCanvas()
    C.cd()
    C.SetLogz()
    H.Draw("col")
    csvfilename = str(sys.argv[1])
    csvfilename = csvfilename.rstrip(".csv")
    C.Print(csvfilename+str(n)+".root")
    # array = H_rotated.GetArray()
    # for y in range(n_bins,0,-1):
    #     # for x in range(17,0,-1):
    #     for x in range(1, (n_bins+1)):
    #         index = x  + y * (n_bins+2)
    #         print(f"{array[index]:.2f}", end=" ")
    #     print("")

import csv
import sys

with open(sys.argv[1]) as csvfile:
    reader = csv.reader(csvfile)
    n = 0
    for row in reader:
        n+=1
        # if 31 < int(float(row[0])) < 41:
        #     print(n)
        if n == int(sys.argv[2]):
            makeplot(row, n)