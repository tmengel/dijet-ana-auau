void count_njets()
{
    // Find one file in the directory
    const char* pattern = "TREE_DIJET_v10_*.root";
    void* dirp = gSystem->OpenDirectory("/sphenix/tg/tg01/jets/dlis/data/v101/runs");
    const char* fname = nullptr;
    while ((fname = gSystem->GetDirEntry(dirp))) {
        if (gSystem->MatchPattern(fname, pattern)) {
            std::cout << "Opening file: " << fname << std::endl;
            break;
        }
    }
    if (!fname) {
        std::cerr << "No file found matching pattern: " << pattern << std::endl;
        return;
    }

    // Open file and get tree
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    TTree* ttree = (TTree*)f->Get("ttree");
    if (!ttree) {
        std::cerr << "No TTree named 'ttree' found!" << std::endl;
        return;
    }

    // Draw histogram of number of jets per event (subtracted jets)
    ttree->Draw("Length$(jet_pt_3_sub)>>hNjets(15,-0.5,14.5)");

    // Style and display
    TH1* hNjets = (TH1*)gDirectory->Get("hNjets");
    if (hNjets) {
        hNjets->SetTitle("Number of jets per event;N_{jets};Events");
        hNjets->SetLineWidth(2);
        hNjets->Draw();
    }
}

