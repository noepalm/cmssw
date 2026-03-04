void plotRecHitDump() {
    
    const char* baselineFile = "mtd_tracking_rechit_histograms.root";
    const char* mergedFile = "mtd_mergedcluster_tracking_rechit_histograms.root";
    
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kCool);
    
    // Output directory
    const char* outDir = "/eos/user/p/pakrap/www/MTD/Clustering/MergedClusters/SinglePi/MTDTrackingRecHit";
    gSystem->mkdir(outDir, kTRUE);
    
    // Open input ROOT files
    TFile* fBaseline = TFile::Open(baselineFile);
    if (!fBaseline || fBaseline->IsZombie()) {
        std::cout << "Error: Cannot open baseline file " << baselineFile << std::endl;
        return;
    }
    
    TFile* fMerged = TFile::Open(mergedFile);
    if (!fMerged || fMerged->IsZombie()) {
        std::cout << "Error: Cannot open merged file " << mergedFile << std::endl;
        return;
    }
    
    std::cout << "Reading histograms from:\n";
    std::cout << "  Baseline: " << baselineFile << "\n";
    std::cout << "  Merged:   " << mergedFile << std::endl;
    
    // Read baseline histograms
    TH1F* h_nHits_base = (TH1F*)fBaseline->Get("mtdTrackingRecHitDump/nHits");
    TH1F* h_time_base = (TH1F*)fBaseline->Get("mtdTrackingRecHitDump/time");
    TH1F* h_timeErr_base = (TH1F*)fBaseline->Get("mtdTrackingRecHitDump/timeErr");
    TH1F* h_energy_base = (TH1F*)fBaseline->Get("mtdTrackingRecHitDump/energy");
    TH1F* h_posX_base = (TH1F*)fBaseline->Get("mtdTrackingRecHitDump/posX");
    TH1F* h_posY_base = (TH1F*)fBaseline->Get("mtdTrackingRecHitDump/posY");
    TH1F* h_errXX_base = (TH1F*)fBaseline->Get("mtdTrackingRecHitDump/errXX");
    TH1F* h_errYY_base = (TH1F*)fBaseline->Get("mtdTrackingRecHitDump/errYY");
    TH2F* h_posXY_base = (TH2F*)fBaseline->Get("mtdTrackingRecHitDump/posXY");
    
    // Read merged histograms
    TH1F* h_nHits_merged = (TH1F*)fMerged->Get("mtdTrackingRecHitDump/nHits");
    TH1F* h_time_merged = (TH1F*)fMerged->Get("mtdTrackingRecHitDump/time");
    TH1F* h_timeErr_merged = (TH1F*)fMerged->Get("mtdTrackingRecHitDump/timeErr");
    TH1F* h_energy_merged = (TH1F*)fMerged->Get("mtdTrackingRecHitDump/energy");
    TH1F* h_posX_merged = (TH1F*)fMerged->Get("mtdTrackingRecHitDump/posX");
    TH1F* h_posY_merged = (TH1F*)fMerged->Get("mtdTrackingRecHitDump/posY");
    TH1F* h_errXX_merged = (TH1F*)fMerged->Get("mtdTrackingRecHitDump/errXX");
    TH1F* h_errYY_merged = (TH1F*)fMerged->Get("mtdTrackingRecHitDump/errYY");
    TH2F* h_posXY_merged = (TH2F*)fMerged->Get("mtdTrackingRecHitDump/posXY");
    
    // Helper function to save canvas
    auto savePlot = [&](TCanvas* c, const char* name) {
        TString pngPath = TString::Format("%s/%s.png", outDir, name);
        TString pdfPath = TString::Format("%s/%s.pdf", outDir, name);
        c->SaveAs(pngPath);
        c->SaveAs(pdfPath);
        std::cout << "Saved " << name << std::endl;
    };
    
    // Plot nHits comparison
    if (h_nHits_base && h_nHits_merged) {
        TCanvas* c = new TCanvas("c_nHits", "Number of Hits", 800, 600);
        gPad->SetLogy();
        h_nHits_base->SetLineColor(kGray+2);
        h_nHits_base->SetLineWidth(2);
        h_nHits_merged->SetLineColor(kMagenta-6);
        h_nHits_merged->SetLineWidth(2);
        
        double max1 = h_nHits_base->GetMaximum();
        double max2 = h_nHits_merged->GetMaximum();
        h_nHits_base->SetMaximum(1.1 * TMath::Max(max1, max2));
        h_nHits_base->GetXaxis()->SetRangeUser(0, 100);
        
        h_nHits_base->Draw("HIST E");
        h_nHits_merged->Draw("HIST E SAME");
        
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(h_nHits_base, "Baseline", "l");
        leg->AddEntry(h_nHits_merged, "Merged Cluster", "l");
        leg->Draw();
        
        savePlot(c, "nHits");
        delete c;
    }
    
    // Plot time comparison
    if (h_time_base && h_time_merged) {
        TCanvas* c = new TCanvas("c_time", "Time", 800, 600);
        h_time_base->SetLineColor(kGray+2);
        h_time_base->SetLineWidth(2);
        h_time_merged->SetLineColor(kMagenta-6);
        h_time_merged->SetLineWidth(2);
        
        double max1 = h_time_base->GetMaximum();
        double max2 = h_time_merged->GetMaximum();
        h_time_base->SetMaximum(1.1 * TMath::Max(max1, max2));
        
        h_time_base->Draw("HIST E");
        h_time_merged->Draw("HIST E SAME");
        
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(h_time_base, "Baseline", "l");
        leg->AddEntry(h_time_merged, "Merged Cluster", "l");
        leg->Draw();
        
        savePlot(c, "time");
        delete c;
    }
    
    // Plot timeErr comparison
    if (h_timeErr_base && h_timeErr_merged) {
        TCanvas* c = new TCanvas("c_timeErr", "Time Error", 800, 600);

        h_timeErr_base->SetLineColor(kGray+2);
        h_timeErr_base->SetLineWidth(2);
        h_timeErr_merged->SetLineColor(kMagenta-6);
        h_timeErr_merged->SetLineWidth(2);
        
        double max1 = h_timeErr_base->GetMaximum();
        double max2 = h_timeErr_merged->GetMaximum();
        h_timeErr_base->SetMaximum(1.1 * TMath::Max(max1, max2));
        
        h_timeErr_base->Draw("HIST E");
        h_timeErr_merged->Draw("HIST E SAME");
        
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(h_timeErr_base, "Baseline", "l");
        leg->AddEntry(h_timeErr_merged, "Merged Cluster", "l");
        leg->Draw();
        
        savePlot(c, "timeErr");
        delete c;
    }
    
    // Plot energy comparison
    if (h_energy_base && h_energy_merged) {
        TCanvas* c = new TCanvas("c_energy", "Energy", 800, 600);
        gPad->SetLogy();

        h_energy_base->SetLineColor(kGray+2);
        h_energy_base->SetLineWidth(2);
        h_energy_merged->SetLineColor(kMagenta-6);
        h_energy_merged->SetLineWidth(2);
        
        double max1 = h_energy_base->GetMaximum();
        double max2 = h_energy_merged->GetMaximum();
        h_energy_base->SetMaximum(1.1 * TMath::Max(max1, max2));
        
        h_energy_base->Draw("HIST E");
        h_energy_merged->Draw("HIST E SAME");
        
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(h_energy_base, "Baseline", "l");
        leg->AddEntry(h_energy_merged, "Merged Cluster", "l");
        leg->Draw();
        
        savePlot(c, "energy");
        delete c;
    }
    
    // Plot posX comparison
    if (h_posX_base && h_posX_merged) {
        TCanvas* c = new TCanvas("c_posX", "Position X", 800, 600);
        h_posX_base->SetLineColor(kGray+2);
        h_posX_base->SetLineWidth(2);
        h_posX_merged->SetLineColor(kMagenta-6);
        h_posX_merged->SetLineWidth(2);
        
        double max1 = h_posX_base->GetMaximum();
        double max2 = h_posX_merged->GetMaximum();
        h_posX_base->SetMaximum(1.1 * TMath::Max(max1, max2));
        h_posX_base->GetXaxis()->SetRangeUser(-6, 6);
        
        h_posX_base->Draw("HIST E");
        h_posX_merged->Draw("HIST E SAME");
        
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(h_posX_base, "Baseline", "l");
        leg->AddEntry(h_posX_merged, "Merged Cluster", "l");
        leg->Draw();
        
        savePlot(c, "posX");
        delete c;
    }
    
    if (h_posY_base && h_posY_merged) {
        TCanvas* c = new TCanvas("c_posY", "Position Y", 800, 600);

        h_posY_base->SetLineColor(kGray+2);
        h_posY_base->SetLineWidth(2);
        h_posY_merged->SetLineColor(kMagenta-6);
        h_posY_merged->SetLineWidth(2);
        
        double max1 = h_posY_base->GetMaximum();
        double max2 = h_posY_merged->GetMaximum();
        h_posY_base->SetMaximum(1.1 * TMath::Max(max1, max2));
        h_posY_base->GetXaxis()->SetRangeUser(-6, 6);
        
        h_posY_base->Draw("HIST E");
        h_posY_merged->Draw("HIST E SAME");
        
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(h_posY_base, "Baseline", "l");
        leg->AddEntry(h_posY_merged, "Merged Cluster", "l");
        leg->Draw();
        
        savePlot(c, "posY");
        delete c;
    }
    
    if (h_errXX_base && h_errXX_merged) {
        TCanvas* c = new TCanvas("c_errXX", "Position Error XX", 800, 600);
        gPad->SetLogy();

        h_errXX_base->SetLineColor(kGray+2);
        h_errXX_base->SetLineWidth(2);
        h_errXX_merged->SetLineColor(kMagenta-6);
        h_errXX_merged->SetLineWidth(2);
        
        double max1 = h_errXX_base->GetMaximum();
        double max2 = h_errXX_merged->GetMaximum();
        h_errXX_base->SetMaximum(1.1 * TMath::Max(max1, max2));
        
        h_errXX_base->Draw("HIST E");
        h_errXX_merged->Draw("HIST E SAME");
        
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(h_errXX_base, "Baseline", "l");
        leg->AddEntry(h_errXX_merged, "Merged Cluster", "l");
        leg->Draw();
        
        savePlot(c, "errXX");
        delete c;
    }
    
    if (h_errYY_base && h_errYY_merged) {
        TCanvas* c = new TCanvas("c_errYY", "Position Error YY", 800, 600);
        gPad->SetLogy();

        h_errYY_base->SetLineColor(kGray+2);
        h_errYY_base->SetLineWidth(2);
        h_errYY_merged->SetLineColor(kMagenta-6);
        h_errYY_merged->SetLineWidth(2);
        
        double max1 = h_errYY_base->GetMaximum();
        double max2 = h_errYY_merged->GetMaximum();
        h_errYY_base->SetMaximum(1.1 * TMath::Max(max1, max2));
        
        h_errYY_base->Draw("HIST E");
        h_errYY_merged->Draw("HIST E SAME");
        
        TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(h_errYY_base, "Baseline", "l");
        leg->AddEntry(h_errYY_merged, "Merged Cluster", "l");
        leg->Draw();
        
        savePlot(c, "errYY");
        delete c;
    }
    if (h_posXY_base) {
        TCanvas* c = new TCanvas("c_posXY_baseline", "Position XY Baseline", 800, 600);
        h_posXY_base->Draw("COLZ");
        savePlot(c, "posXY_baseline");
        delete c;
    }
    
    if (h_posXY_merged) {
        TCanvas* c = new TCanvas("c_posXY_merged", "Position XY Merged", 800, 600);
        h_posXY_merged->Draw("COLZ");
        savePlot(c, "posXY_merged");
        delete c;
    }
    
    fBaseline->Close();
    fMerged->Close();
    
    std::cout << "\nAll plots saved to: " << outDir << std::endl;
}