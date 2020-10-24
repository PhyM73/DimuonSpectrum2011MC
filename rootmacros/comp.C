#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"

void comp(){
TFile *f=new TFile("DM2011_data.root"); // Data File
TH1F *h= (TH1F*)f->Get("demo/GM_mass_tight_iso");

/***** Rescale *****/
double lumi = 2.11; // inverse femtobarn

TFile *f1=new TFile("DMDYMC_accept.root"); // DMDYMC File
TFile *f2=new TFile("DMZMC_accept.root"); // ZMC  File

TH1F *h0= (TH1F*)f1->Get("demo/Mmultiplicity");
TH1F *hz0= (TH1F*)f2->Get("demo/Mmultiplicity");

double adjustdy = 1.2;
double adjustz  = 1.2;
double xsecdy = 9507000. ; // femtobarn
double xsecz  = 2475000. ; // femtobarn 
double rdy = xsecdy*lumi/(h0->GetEntries())*adjustdy;
double rz = xsecz*lumi/(hz0->GetEntries())*adjustz;

cout<<setprecision(5)<<"rescale DYMC: "<<rdy<<endl;
cout<<setprecision(5)<<"rescale ZMC:  "<<rz<<endl;

TH1F *h1= (TH1F*)f1->Get("demo/GM_mass_tight_iso");
TH1F *hz1= (TH1F*)f2->Get("demo/GM_mass_tight_iso");

h1->Scale(rdy);
hz1->Scale(rz);
//h1->Scale(rdy,"nosw2");
//hz1->Scale(rz,"nosw2");

TH1F *he = (TH1F*) hz1->Clone();
hz1->Add(h1);
he->Add(h1);


/***** Divide *****/
TH1F *hr1 = (TH1F*) hz1-> Clone();
//hr1->Sumw2();
hr1->Divide(hz1);
TH1F *hre = (TH1F*) hr1->Clone();

TH1F *hratio = (TH1F*) h-> Clone();
//hratio->Sumw2();
hratio->Divide(hz1);


/***** Draw *****/
TCanvas *c= new TCanvas("c","GM_mass",580,725);
gStyle->SetOptStat("");

TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
pad1->SetBottomMargin(0.15); 
pad1->Draw();  // Draw the upper pad: pad1
pad1->cd();    // pad1 becomes the current pad

he->Draw("E2,][");
hz1->Draw("HIST,SAME,][");
h->Draw("E1,SAME,][");

hz1->SetLineColor(kOrange+7);
hz1->SetMarkerColor(kOrange+7);
he->SetFillColor(kOrange-9);
he->SetMarkerColor(kOrange+7);
he->SetLineColor(kOrange+7);
h->SetLineColor(kBlack);

he->SetTitle("");
he->SetXTitle("m_{#mu#mu}[GeV]");
he->GetXaxis()->SetTitleSize(18);
he->GetXaxis()->SetTitleFont(45);
he->GetXaxis()->SetTitleOffset(1.6);
he->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
he->GetXaxis()->SetLabelSize(15);
he->GetXaxis()->SetRangeUser(16,150);

he->SetYTitle("");
he->SetMinimum(10.);
he->SetMaximum(310000.);
he->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
he->GetYaxis()->SetLabelSize(15);
pad1->SetLogy();


TLegend* l = new TLegend(0.9,0.9,0.7,0.7);
l->AddEntry(h,"CMS11a","lep");
l->AddEntry(he,"CMS MC");
l->Draw("SAME");

c->cd();          // Go back to the main canvas before defining pad2

TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->SetTopMargin(0.05);
pad2->SetBottomMargin(0.3);
pad2->Draw();
pad2->cd();       // pad2 becomes the current pad


hre->Draw("E2 ][");
hr1->Draw("HIST SAME ][");
hratio->Draw("E1 SAME ][");

hr1->SetLineColor(kOrange+5);
hr1->SetLineWidth(2);
hre->SetFillColor(kOrange-9);
hre->SetMarkerColor(kOrange-9);
hratio->SetLineColor(kBlack);

hre->SetMinimum(0.5);
hre->SetMaximum(1.5);
hre->SetTitle("");
hre->SetXTitle("m_{#mu#mu}[GeV]");
hre->GetXaxis()->SetTitleSize(18);
hre->GetXaxis()->SetTitleFont(45);
hre->GetXaxis()->SetTitleOffset(4);
hre->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
hre->GetXaxis()->SetLabelSize(15);
hre->GetXaxis()->SetRangeUser(16,150);

hre->SetYTitle("");
hre->SetNdivisions(505,"y");
hre->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
hre->GetYaxis()->SetLabelSize(15);

}
