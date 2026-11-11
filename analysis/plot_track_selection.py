#!/usr/bin/env python3
"""Plot useful yield and provisional Ay-squared information by track policy."""
import argparse
import csv
import json
from array import array
from pathlib import Path
import ROOT


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('comparison_stem')
    args=parser.parse_args()
    result=json.loads(Path(args.comparison_stem+'.json').read_text())
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(.15)
    ROOT.gStyle.SetPadBottomMargin(.13)
    policies=('legacy','single','min_angle_uncut','min_angle')
    labels=('Original leading-momentum reference','Exactly one front-pair track (no 1 deg cut)',
            'Previous smallest-angle rule (no 1 deg cut)',
            f'Smallest-angle track with angle > {result["min_angle_strict_lower_bound_deg"]:g} deg')
    colors=(ROOT.kGray+2,ROOT.kBlue+1,ROOT.kGreen+2,ROOT.kRed+1)
    hist={}
    for p in policies:
        hist[p]=[ROOT.TH1D(p+'_yield',';Polar scattering angle (deg);Useful events / incident proton / 0.25 deg',40,0,10),
                 ROOT.TH1D(p+'_fom',';Polar scattering angle (deg);Intrinsic FOM / incident proton / 0.25 deg',40,0,10)]
        for h in hist[p]: h.Sumw2()
    with open(args.comparison_stem+'_events.csv') as f:
        for row in csv.DictReader(f):
            for p in policies:
                if int(row[p+'_useful']):
                    angle=float(row[p+'_theta_deg']); weight=float(row['weight'])
                    hist[p][0].Fill(angle,weight)
                    hist[p][1].Fill(angle,weight*float(row[p+'_Ay2']))
    canvas=ROOT.TCanvas('track_policies','FPP track policies',1500,750)
    canvas.Divide(2,1)
    legends=[]
    titles=[]
    for panel in range(2):
        canvas.cd(panel+1)
        maximum=max(hist[p][panel].GetMaximum() for p in policies)/result['incident_weight']
        legend=ROOT.TLegend(.18,.67,.94,.89)
        legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(.025)
        for index,p in enumerate(policies):
            h=hist[p][panel]; h.Scale(1./result['incident_weight'])
            h.SetLineColor(colors[index]); h.SetLineWidth(3)
            h.GetYaxis().SetTitleOffset(1.7)
            h.SetMaximum(maximum*1.42); h.SetMinimum(0)
            h.Draw('hist' if index==0 else 'hist same')
            legend.AddEntry(h,labels[index],'l')
        legend.Draw(); legends.append(legend)
        title=ROOT.TLatex(); title.SetNDC(); title.SetTextSize(.03)
        title.DrawLatex(.16,.94,'Useful-p_{T} yield' if panel==0 else 'A_{y}^{2}-weighted information (multi-track provisional)')
        titles.append(title)
    canvas.SaveAs(args.comparison_stem+'.png')
    canvas.SaveAs(args.comparison_stem+'.pdf')


if __name__=='__main__': main()
