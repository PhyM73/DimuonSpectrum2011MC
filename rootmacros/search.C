#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"

void search(){
TFile *f=new TFile("DM2011_data.root"); // Data File
TH1F *h11= (TH1F*)f->Get("demo/GM_mass_iso_0");
TH1F *h12= (TH1F*)f->Get("demo/GM_mass_iso_25");
TH1F *h13= (TH1F*)f->Get("demo/GM_mass_iso_60");
TH1F *h21= (TH1F*)f->Get("demo/GM_mass_pro_0");
TH1F *h22= (TH1F*)f->Get("demo/GM_mass_pro_25");
TH1F *h23= (TH1F*)f->Get("demo/GM_mass_pro_60");


/***** Draw *****/
TCanvas *c= new TCanvas("c","GM_mass",1130,538);
gStyle->SetOptStat("");

TPad *pad1 = new TPad("pad1","pad1",0,0,0.5,1.0);
pad1->SetBottomMargin(0.2); 
pad1->Draw();  // Draw the upper pad: pad1
pad1->cd();    // pad1 becomes the current pad

h11->Draw("E1,][");
h12->Draw("E1,SAME,][");
h13->Draw("E1,SAME,][");

h11->SetLineColor(kBlue);
h11->SetMarkerColor(kBlue);
h12->SetLineColor(kBlack);
h12->SetMarkerColor(kBlack);
h13->SetLineColor(kGreen);
h13->SetMarkerColor(kGreen);

h11->SetTitle("");
h11->SetXTitle("m_{#mu#mu}[GeV]");
h11->GetXaxis()->SetTitleSize(18);
h11->GetXaxis()->SetTitleFont(45);
h11->GetXaxis()->SetTitleOffset(1.6);
h11->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
h11->GetXaxis()->SetLabelSize(15);

h11->SetYTitle("");
h11->SetMinimum(80.);
h11->SetMaximum(200000.);
h11->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
h11->GetYaxis()->SetLabelSize(15);
pad1->SetLogy();


TLegend* l1 = new TLegend(0.9,0.9,0.6,0.7);
l1->AddEntry(h11,"no p^{#mu#mu}_{T} cut","lep");
l1->AddEntry(h12,"p^{#mu#mu}_{T}>25Gev/c","lep");
l1->AddEntry(h13,"p^{#mu#mu}_{T}>60Gev/c","lep");
l1->Draw("SAME");

c->cd();          // Go back to the main canvas before defining pad2

TPad *pad2 = new TPad("pad2", "pad2", 0.5, 0, 1., 1.);
//pad2->SetLeftMargin(0.05);
pad2->SetBottomMargin(0.2);
pad2->Draw();
pad2->cd();       // pad2 becomes the current pad

h21->Draw("E1,][");
h22->Draw("E1,SAME,][");
h23->Draw("E1,SAME,][");

h21->SetLineColor(kBlue);
h21->SetMarkerColor(kBlue);
h22->SetLineColor(kBlack);
h22->SetMarkerColor(kBlack);
h23->SetLineColor(kGreen);
h23->SetMarkerColor(kGreen);

h21->SetTitle("");
h21->SetXTitle("m_{#mu#mu}[GeV]");
h21->GetXaxis()->SetTitleSize(18);
h21->GetXaxis()->SetTitleFont(45);
h21->GetXaxis()->SetTitleOffset(1.6);
h21->GetXaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
h21->GetXaxis()->SetLabelSize(15);

h21->SetYTitle("");
h21->SetMinimum(80.);
h21->SetMaximum(200000.);
h21->GetYaxis()->SetLabelFont(45); // Absolute font size in pixel (precision 3)
h21->GetYaxis()->SetLabelSize(15);
pad2->SetLogy();


TLegend* l2 = new TLegend(0.9,0.9,0.6,0.7);
l2->AddEntry(h21,"no p^{#mu#mu}_{T} cut","lep");
l2->AddEntry(h22,"p^{#mu#mu}_{T}>25Gev/c","lep");
l2->AddEntry(h23,"p^{#mu#mu}_{T}>60Gev/c","lep");
l2->Draw("SAME");

}
