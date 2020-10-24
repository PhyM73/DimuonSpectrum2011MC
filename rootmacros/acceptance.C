#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"

void acceptance(){

TFile *f=new TFile("DMZMC_accept.root"); // Data File
TH1F *h= (TH1F*)f->Get("demo/Accept");

for (int i=1; i<3; ++i){
cout<<fixed<<setprecision(0)<<h->GetBinContent(i)<<endl;
}
cout<<setprecision(4)<<(h->GetBinContent(2))/(h->GetBinContent(1))<<endl;
}
