#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TPaveStats.h>

/*
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CUDADataFormats/HcalDigi/interface/DigiCollection.h"
*/

#define CREATE_HIST_1D(varname, nbins, first, last) \
    auto varname = new TH1D(#varname, #varname, nbins, first, last)

#define CREATE_HIST_2D(varname, nbins, first, last) \
    auto varname = new TH2D(#varname, #varname, nbins, first, last, nbins, first, last)

int main(int argc, char *argv[]) {
    if (argc<3) {
        std::cout << "run with: ./<exe> <path to input file> <path to output file>\n";
        exit(0);
    }

    std::string inFileName{argv[1]};
    std::string outFileName{argv[2]};
    
    // branches to use
    /*
    edm::Wrapper<edm::DataFrameContainer> *wgpuf01he=nullptr, *wgpuf5hb=nullptr;
    edm::Wrapper<QIE11DigiCollection> *wcpuf01he=nullptr;
    edm::Wrapper<HBHEDigiCollection> *wcpuf5hb=nullptr;

    std::string inFileName{argv[1]};
    std::string outFileName{argv[2]};
    */

    // prep output 
    TFile rfout{outFileName.c_str(), "recreate"};

    /*
    CREATE_HIST_1D(hADCf01HEGPU, 256, 0, 256);
    CREATE_HIST_1D(hADCf01HECPU, 256, 0, 256);
    CREATE_HIST_1D(hADCf5HBGPU, 128, 0, 128);
    CREATE_HIST_1D(hADCf5HBCPU, 128, 0, 128);
    CREATE_HIST_1D(hTDCf01HEGPU, 64, 0, 64);
    CREATE_HIST_1D(hTDCf01HECPU, 64, 0, 64);

    CREATE_HIST_2D(hADCf01HEGPUvsCPU, 256, 0, 256);
    CREATE_HIST_2D(hADCf5HBGPUvsCPU, 128, 0, 128);
    CREATE_HIST_2D(hTDCf01HEGPUvsCPU, 64, 0, 64);
    */

    /*
    auto hADCEBGPU = new TH1D("hADCEBGPU", "hADCEBGPU", nbins, 0, last);
    auto hADCEBCPU = new TH1D("hADCEBCPU", "hADCEBCPU", nbins, 0, last);
    auto hADCEEGPU = new TH1D("hADCEEGPU", "hADCEEGPU", nbins, 0, last);
    auto hADCEECPU = new TH1D("hADCEECPU", "hADCEECPU", nbins, 0, last);

    auto hGainEBGPU = new TH1D("hGainEBGPU", "hGainEBGPU", 4, 0, 4);
    auto hGainEBCPU = new TH1D("hGainEBCPU", "hGainEBCPU", 4, 0, 4);
    auto hGainEEGPU = new TH1D("hGainEEGPU", "hGainEEGPU", 4, 0, 4);
    auto hGainEECPU = new TH1D("hGainEECPU", "hGainEECPU", 4, 0, 4);

    auto hADCEBGPUvsCPU = new TH2D("hADCEBGPUvsCPU", "hADCEBGPUvsCPU",
        nbins, 0, last, nbins, 0, last);
    auto hADCEEGPUvsCPU = new TH2D("hADCEEGPUvsCPU", "hADCEEGPUvsCPU",
        nbins, 0, last, nbins, 0, last);
    auto hGainEBGPUvsCPU = new TH2D("hGainEBGPUvsCPU", "hGainEBGPUvsCPU",
        4, 0, 4, 4, 0, 4);
    auto hGainEEGPUvsCPU = new TH2D("hGainEEGPUvsCPU", "hGainEEGPUvsCPU",
        4, 0, 4, 4, 0, 4);
        */

    // prep input
    TFile rfin{inFileName.c_str()};
    TTree *rt = (TTree*)rfin.Get("Events");
    /*
    rt->SetBranchAddress("QIE11DataFrameHcalDataFrameContainer_hcalDigis__RECO.",
        &wcpuf01he);
    rt->SetBranchAddress("HBHEDataFramesSorted_hcalDigis__RECO.", &wcpuf5hb);
    rt->SetBranchAddress("edmDataFrameContainer_hcalCPUDigisProducer_f5HBDigis_RECO.",
        &wgpuf5hb);
    rt->SetBranchAddress(
        "edmDataFrameContainer_hcalCPUDigisProducer_f01HEDigis_RECO.",
        &wgpuf01he);
        */

    // accumulate
    auto const nentries = rt->GetEntries();
    std::cout << ">>> nentries = " << nentries << std::endl;
    for (int ie=0; ie<nentries; ++ie) {
        rt->GetEntry(ie);

        /*
        auto const& f01HEProduct = wgpuf01he->bareProduct();
        auto const& f5HBProduct = wgpuf5hb->bareProduct();
        auto const& qie11Product = wcpuf01he->bareProduct();
        auto const qie11Filtered = filterQIE11(qie11Product);
        auto const& qie8Product = wcpuf5hb->bareProduct();

        auto const ngpuf01he = f01HEProduct.size();
        auto const ngpuf5hb = f5HBProduct.size();
        auto const ncpuf01he = qie11Product.size();
        auto const ncpuf5hb = qie8Product.size();

        printf("ngpuf01he = %u nqie11 = %u ncpuf01he = %u ngpuf5hb = %u ncpuf5hb = %u\n",
            f01HEProduct.size(), qie11Product.size(), qie11Filtered.size(), 
            f5HBProduct.size(),
            static_cast<uint32_t>(qie8Product.size()));

        if (ngpuf01he != ncpuf01he) {
            std::cerr << "*** mismatch in number of flavor 01 digis for event "
                      << ie
                      << std::endl
                      << ">>> ngpuf01he = " << ngpuf01he << std::endl
                      << ">>> ncpuf01he = " << ncpuf01he
                      << std::endl;
        }

        {
            auto const& idsgpu = f01HEProduct.ids();
            auto const& datagpu = f01HEProduct.data();

            for (uint32_t ich=0; ich<ncpuf01he; ich++) {
                auto const cpudf = QIE11DataFrame{qie11Product[ich]};
                auto const cpuid = cpudf.id();
                auto iter2idgpu = std::find(idsgpu.begin(), idsgpu.end(), cpuid);

                if (iter2idgpu == idsgpu.end()) {
                    std::cerr << "missing "<< HcalDetId{cpuid} << std::endl;
                    continue;
                }

                // FIXME: cna fail...
                assert(*iter2idgpu == cpuid);

                auto const ptrdiff = iter2idgpu - idsgpu.begin();
                auto const nsamples_gpu = hcal::compute_nsamples<hcal::Flavor01>(
                    f01HEProduct.stride());
                auto const nsamples_cpu = qie11Product.samples();
                assert(static_cast<uint32_t>(nsamples_cpu) == nsamples_gpu);

                uint32_t ichgpu = ptrdiff;
                uint32_t offset = ichgpu * f01HEProduct.stride();
                uint16_t const* df_start = datagpu.data() + offset;
                for (uint32_t sample=0u; sample<nsamples_gpu; sample++) {
                    auto const cpuadc = cpudf[sample].adc();
                    auto const gpuadc = hcal::adc_for_sample<hcal::Flavor01>(
                        df_start, sample);
                    auto const cputdc = cpudf[sample].tdc();
                    auto const gputdc = hcal::tdc_for_sample<hcal::Flavor01>(
                        df_start, sample);

                    hADCf01HEGPU->Fill(gpuadc);
                    hADCf01HECPU->Fill(cpuadc);
                    hTDCf01HEGPU->Fill(gputdc);
                    hTDCf01HECPU->Fill(cputdc);
                    hADCf01HEGPUvsCPU->Fill(cpuadc, gpuadc);
                    hTDCf01HEGPUvsCPU->Fill(cputdc, gputdc);

                    // At RAW Decoding level there must not be any mistmatches 
                    // in the adc values at all!
                    assert(static_cast<uint8_t>(cpuadc) == gpuadc);
                    assert(static_cast<uint8_t>(cputdc) == gputdc);
                }
            }
        }

        if (ngpuf5hb != ncpuf5hb) {
            std::cerr << "*** mismatch in number of flavor 5 digis for event "
                      << ie
                      << std::endl
                      << ">>> ngpuf5hb = " << ngpuf5hb << std::endl
                      << ">>> ncpuf5hb = " << ncpuf5hb
                      << std::endl;
        }

        {
            auto const& idsgpu = f5HBProduct.ids();
            auto const& datagpu = f5HBProduct.data();
            for (uint32_t i=0; i<ncpuf5hb; i++) {
                auto const cpudf = qie8Product[i];
                auto const cpuid = cpudf.id().rawId();
                auto iter2idgpu = std::find(idsgpu.begin(), idsgpu.end(), cpuid);
                if (iter2idgpu == idsgpu.end()) {
                    std::cerr << "missing "<< HcalDetId{cpuid} << std::endl;
                    continue;
                }

                assert(*iter2idgpu == cpuid);

                auto const ptrdiff = iter2idgpu - idsgpu.begin();
                auto const nsamples_gpu = hcal::compute_nsamples<hcal::Flavor5>(
                    f5HBProduct.stride());
                auto const nsamples_cpu = qie8Product[0].size();
                assert(static_cast<uint32_t>(nsamples_cpu) == nsamples_gpu);

                uint32_t offset = ptrdiff * f5HBProduct.stride();
                uint16_t const* df_start = datagpu.data() + offset;
                for (uint32_t sample=0u; sample<nsamples_gpu; sample++) {
                    auto const cpuadc = cpudf.sample(sample).adc();
                    auto const gpuadc = hcal::adc_for_sample<hcal::Flavor5>(
                        df_start, sample);

                    hADCf5HBGPU->Fill(gpuadc);
                    hADCf5HBCPU->Fill(cpuadc);
                    hADCf5HBGPUvsCPU->Fill(cpuadc, gpuadc);

                    // the must for us at RAW Decoding stage
                    assert(static_cast<hcal::Flavor5::adc_type>(cpuadc) == gpuadc);
                }
            }
        }
        */
    }

    /*
    {
        TCanvas c{"plots", "plots", 4200, 6200};
        c.Divide(2, 3);
        c.cd(1);
        {
            gPad->SetLogy();
            hADCf01HECPU->SetLineColor(kBlack);
            hADCf01HECPU->SetLineWidth(1.);
            hADCf01HECPU->Draw("");
            hADCf01HEGPU->SetLineColor(kBlue);
            hADCf01HEGPU->SetLineWidth(1.);
            hADCf01HEGPU->Draw("sames");
            gPad->Update();
            auto stats = (TPaveStats*)hADCf01HEGPU->FindObject("stats");
            auto y2 = stats->GetY2NDC();
            auto y1 = stats->GetY1NDC();
            stats->SetY2NDC(y1);
            stats->SetY1NDC(y1 - (y2-y1));
        }
        c.cd(2);
        {
            gPad->SetLogy();
            hADCf5HBCPU->SetLineColor(kBlack);
            hADCf5HBCPU->SetLineWidth(1.);
            hADCf5HBCPU->Draw("");
            hADCf5HBGPU->SetLineColor(kBlue);
            hADCf5HBGPU->SetLineWidth(1.);
            hADCf5HBGPU->Draw("sames");
            gPad->Update();
            auto stats = (TPaveStats*)hADCf5HBGPU->FindObject("stats");
            auto y2 = stats->GetY2NDC();
            auto y1 = stats->GetY1NDC();
            stats->SetY2NDC(y1);
            stats->SetY1NDC(y1 - (y2-y1));
        }
        c.cd(3);
        hADCf01HEGPUvsCPU->Draw("colz");
        c.cd(4);
        hADCf5HBGPUvsCPU->Draw("colz");
        c.cd(5);
        {
            gPad->SetLogy();
            hTDCf01HECPU->SetLineColor(kBlack);
            hTDCf01HECPU->SetLineWidth(1.);
            hTDCf01HECPU->Draw("");
            hTDCf01HEGPU->SetLineColor(kBlue);
            hTDCf01HEGPU->SetLineWidth(1.);
            hTDCf01HEGPU->Draw("sames");
            gPad->Update();
            auto stats = (TPaveStats*)hTDCf01HEGPU->FindObject("stats");
            auto y2 = stats->GetY2NDC();
            auto y1 = stats->GetY1NDC();
            stats->SetY2NDC(y1);
            stats->SetY1NDC(y1 - (y2-y1));
        }
        c.cd(6);
        hTDCf01HEGPUvsCPU->Draw("colz");
        
        c.SaveAs("plots.pdf");
    }*/

    rfin.Close();
    rfout.Write();
    rfout.Close();
}
