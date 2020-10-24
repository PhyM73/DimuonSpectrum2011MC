#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"

void cutflow(){

TFile *f=new TFile("DM2011_data.root"); // Data File
TH1F *h= (TH1F*)f->Get("demo/Cut_Flow");

for (int i=1; i<15; ++i){
cout<<fixed<<setprecision(0)<<h->GetBinContent(i)<<endl;
}

}
