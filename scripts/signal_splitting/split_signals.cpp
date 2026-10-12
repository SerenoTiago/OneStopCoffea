#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>
#include <format>
#include <memory>
#include <fstream>
#include <filesystem>
#include <stdexcept>
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

bool matchesRequestedModel(const std::string& branch_name, const std::unordered_set<std::string>& requested_models){
    if(requested_models.empty()) {return true;}
    return requested_models.contains(branch_name)
        || (branch_name.starts_with("GenModel_") && requested_models.contains(branch_name.substr(9)));
}


void processOneFile(const std::filesystem::path& infile, const std::filesystem::path& outdir, const std::unordered_set<std::string>& requested_models = {}){
    std::cout << std::format("Processing {} file to outdirectory {}\n", infile.string(), outdir.string());

    auto base = infile.stem().string();

    std::vector<std::string> gen_names;
    std::vector<std::string> other_names;

    ROOT::RDataFrame rdf("Events", infile.string());

    for(const auto name: rdf.GetColumnNames()){
        if(name.find("GenModel") == std::string::npos){
            other_names.push_back(name);
        } else {
            if(matchesRequestedModel(name, requested_models)){
                gen_names.push_back(name);
            }

        }
    }

    if(gen_names.empty()){
        throw std::runtime_error(requested_models.empty() ? "No GenModel columns found" : "Could not find any requested GenModel columns");
    }

    for(const auto& name : gen_names){
        auto filtered = rdf.Filter([](bool x){return x;}, {name});
        // Group T2tt source-file shards by the generator-level stop mass.
        // Other models retain the historical layout because their second
        // model-name field is not necessarily a stop mass.
        const auto model_name = name.starts_with("GenModel_") ? name.substr(9) : name;
        auto relative_name = std::filesystem::path(base) / (std::string("signal") + name.substr(8,name.size()-1) + ".root");
        if(model_name.starts_with("T2tt_")){
            const auto first_separator = model_name.find('_');
            const auto second_separator = model_name.find('_', first_separator + 1);
            if(second_separator == std::string::npos){
                throw std::runtime_error(std::format("Could not parse T2tt stop mass from {}", name));
            }
            const auto stop_mass = model_name.substr(first_separator + 1, second_separator - first_separator - 1);
            relative_name = std::filesystem::path("mStop-" + stop_mass) / relative_name;
        }
        const auto final_name = getUniquePath(outdir / relative_name);
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
    if(argc < 3){
        std::cerr << std::format("Usage: {} INPUT OUTPUT_DIR [MODEL ...]\n", argv[0]);
        return 1;
    }
    std::string fname(argv[1]);
    std::string outdir(argv[2]);
    std::unordered_set<std::string> models;
    for(int i = 3; i < argc; ++i){
        models.insert(argv[i]);
    }
    processOneFile(fname, outdir, models);
    return 0;
}
