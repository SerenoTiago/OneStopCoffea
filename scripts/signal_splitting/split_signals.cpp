#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>
#include <format>
#include <memory>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <optional>
#include "TFile.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RSnapshotOptions.hxx"


std::filesystem::path getUniquePath(const std::filesystem::path& p){
    if(! std::filesystem::exists(p)) {return p;}
    auto ret = p;
    auto parent = p.parent_path();
    auto stem = p.stem();
    auto ext = p.extension();
    int i = 1;
    while(std::filesystem::exists(ret)){
        ret  = parent / (stem.string() + "_" + std::to_string(i)) ;
        ret.replace_extension(ext);
        ++i;
    }
    return ret;
}


std::vector<std::string> getFiles(){
    std::vector<std::string> file_paths;
    std::ifstream input("signal_files_2018.txt");
    std::string line;

    while(std::getline(input, line)){
        file_paths.push_back(std::string("root://cmsxrootd.fnal.gov//") + line);
    }
    return file_paths;

}

bool matchesRequestedModel(const std::string& branch_name, const std::optional<std::string>& requested_model){
    if(!requested_model) {return true;}
    return branch_name == *requested_model || branch_name == "GenModel_" + *requested_model;
}


void processOneFile(const std::filesystem::path& infile, const std::filesystem::path& outdir, const std::optional<std::string>& requested_model = std::nullopt){
    std::cout << std::format("Processing {} file to outdirectory {}\n", infile.string(), outdir.string());

    auto base = infile.stem().string();

    std::vector<std::string> gen_names;
    std::vector<std::string> other_names;

    ROOT::RDataFrame rdf("Events", infile.string());

    for(const auto name: rdf.GetColumnNames()){
        if(name.find("GenModel") == std::string::npos){
            other_names.push_back(name);
        } else {
            if(matchesRequestedModel(name, requested_model)){
                gen_names.push_back(name);
            }

        }
    }

    if(gen_names.empty()){
        throw std::runtime_error(requested_model ? std::format("Could not find requested model {}", *requested_model) : "No GenModel columns found");
    }

    for(const auto& name : gen_names){
        auto filtered = rdf.Filter([](bool x){return x;}, {name});
        const auto final_name = getUniquePath(outdir / base / ( std::string("signal") + name.substr(8,name.size()-1) + ".root"));
        std::filesystem::create_directories(final_name.parent_path());
        std::cout << std::format("Saving {} to {}\n", name, final_name.string());
        ROOT::RDF::RSnapshotOptions opts;
        opts.fLazy = false;
        opts.fAutoFlush = 1000;
        opts.fCompressionLevel = 0;
        try {
            filtered.Snapshot("Events", final_name.string(), other_names, opts);
        } catch(const std::exception& e) {
            std::cerr << std::format("Failed while saving {}: {}\n", name, e.what());
            throw;
        }
    }
}


int main(int argc, char* argv[]) {
    if(argc != 3 && argc != 4){
        std::cerr << std::format("Usage: {} INPUT OUTPUT_DIR [MODEL]\n", argv[0]);
        return 1;
    }
    std::string fname(argv[1]);
    std::string outdir(argv[2]);
    std::optional<std::string> model;
    if(argc == 4){
        model = argv[3];
    }
    processOneFile(fname, outdir, model);
    return 0;
}
