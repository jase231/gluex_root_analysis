#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <chrono>
#include <iostream>

#include "compare_hypotheses.h"
#include "inih/INIReader.h"

using namespace std::chrono;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Please pass the configuration file.\nUsage: " << argv[0]
              << " <config.ini>\n";
    return 1;
  }
  // init benchmarking
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  // read in config file
  INIReader reader(argv[1]);
  if (reader.ParseError() != 0) {
    std::cerr << "Error reading config file " << argv[1] << '\n';
    return 1;
  }

  // optional configs
  std::string out_file = reader.Get("Misc", "outfile", "placeholder");
  bool best_by_beam = reader.GetBoolean("Misc", "best_per_beam", true);
  bool preserve_combos = reader.GetBoolean("Misc", "preserve_combos", true);
  bool logging = reader.GetBoolean("Misc", "logging", false);

  std::string tree1 = reader.Get("1", "tree", "");
  std::string glob1 = reader.Get("1", "glob", "");
  if (glob1.empty() || tree1.empty()) {
    std::cerr << "Primary hypothesis parameters missing. Please enter the "
                 "primary hypotheses' filename and treename in the config.\n";
    return 1;
  }
  Tree_config primary = {tree1, glob1};

  // get number of alternative hypotheses
  int num_alt_hypos = reader.GetInteger("1", "num_alt_hypos", 1);
  std::vector<Tree_config> alt_hypo_configs;
  alt_hypo_configs.reserve(num_alt_hypos);
  for (int i = 0; i < num_alt_hypos; i++) {
    std::string num_as_string = std::to_string(i + 2);
    std::string glob = reader.Get(num_as_string, "glob", "");
    std::string tree = reader.Get(num_as_string, "tree", "");

    if (glob.empty() || tree.empty()) {
      std::cerr
          << "At least one alternative hypothesis parameter is missing. Please "
             "ensure you entered the filenames' glob and treenames for the "
             "number of hypotheses you indicated in the config.\n";
      return 1;
    }

    alt_hypo_configs.push_back({glob, tree});
  }

  if (best_by_beam) {
    std::cout << "Running in best combo per beam ID mode.\n";
  } else {
    std::cout << "Running in best overall combo mode.\n";
  }

  std::cout << "Comparing:\n";
  for (Tree_config tree : alt_hypo_configs) {
    std::cout << tree1 + " with " + tree.treename + '\n';
  }

  std::cout << "Pre-processing data..." << std::endl;
  compare_hypotheses c(glob1, tree1, alt_hypo_configs, best_by_beam);
  c.set_preserving(preserve_combos);
  c.set_logging(logging);
  c.set_match_by_beam(best_by_beam);

  c.prepare_data();
  std::cout << "Data prepared, finding matches..." << std::endl;
  c.find_matches();

  std::cout << "Number of matches: " << c.matches << std::endl;

  // benchmark the matching process
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto matching_duration = duration_cast<microseconds>(t2 - t1).count();
  std::cout << "The matching process took: " << matching_duration * 1E-6
            << " seconds\n";

  std::cout << "Writing to file...\n";
  c.write_to_file(out_file);

  // benchmark the writing-to-file
  high_resolution_clock::time_point t3 = high_resolution_clock::now();
  auto writeDuration = duration_cast<microseconds>(t3 - t2).count();
  std::cout << "The file-writing process took: " << writeDuration * 1E-6
            << " seconds\n";

  return 0;
}
