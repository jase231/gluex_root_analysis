#include "compare_hypotheses.h"

#include <cmath>
#include <cstddef>
#include <fstream>

// constructor for the hypothesis trees. passes input directly to RDataFrame
// constructor.
hypothesis_tree_base::hypothesis_tree_base(std::string glob,
                                           std::string tree_name,
                                           bool match_type)
    : df(ROOT::RDataFrame(tree_name, glob)),
      tree_name(tree_name),
      match_by_best_per_beam(match_type),
      logging(false) {}

// same as containsEventID() but uses std::pair as key. to be replaced with
// template.
bool hypothesis_tree_base::contains_event_id(
    std::pair<unsigned long long, unsigned> pair_key) const {
  if (!match_by_best_per_beam) {
    return (event_as_key_map.find(pair_key.first) != event_as_key_map.end());
  }
  return (event_beam_as_key_map.find(pair_key) != event_beam_as_key_map.end());
}

// loads relevant columns from the RDataFrame
void hypothesis_tree_base::fill_column_vecs() {
  event_column_data = *df.Take<unsigned long long>("event");
  run_column_data = *df.Take<unsigned int>("run");
  beam_column_data = *df.Take<unsigned int>("beam_beamid");
  chi_sq_column_data = *df.Take<float>("kin_chisq");
  ndf_column_data = *df.Take<unsigned>("kin_ndf");
}

// fill each combo with data from the data columns
void hypothesis_tree_best_combo::update_combo_data(size_t index) {
  // set eventAsKey map using event ID as key
  event_as_key_map[event_column_data[index]].set_chi_sq(
      chi_sq_column_data[index]);
  event_as_key_map[event_column_data[index]].set_event(
      event_column_data[index]);
  event_as_key_map[event_column_data[index]].set_run(run_column_data[index]);
  event_as_key_map[event_column_data[index]].set_beam(beam_column_data[index]);
  event_as_key_map[event_column_data[index]].set_ndf(ndf_column_data[index]);
}

// set combo data for best combo per beam
void hypothesis_tree_best_per_beam::update_combo_data(size_t index) {
  // set eventBeamAsKey map using event ID and beam ID as key
  auto pair_key =
      std::make_pair(event_column_data[index], beam_column_data[index]);

  event_beam_as_key_map[pair_key].set_chi_sq(chi_sq_column_data[index]);
  event_beam_as_key_map[pair_key].set_event(event_column_data[index]);
  event_beam_as_key_map[pair_key].set_run(run_column_data[index]);
  event_beam_as_key_map[pair_key].set_beam(beam_column_data[index]);
  event_beam_as_key_map[pair_key].set_ndf(ndf_column_data[index]);
  return;
}

// for every row, remove duplicate rows (by shared event ID) and keep only combo
// lowest chisq
void hypothesis_tree_best_combo::filter_high_chi_sq_events() {
  for (size_t i = 0; i < event_column_data.size(); ++i) {
    auto pair_key =
        std::make_pair(event_column_data[i], 1851);  // placeholder beam ID
    if (!contains_event_id(pair_key)) {
      update_combo_data(i);
    } else {
      if (event_as_key_map[event_column_data[i]].get_chi_sq() >
          chi_sq_column_data[i]) {
        update_combo_data(i);
      }
    }
  }
}

void hypothesis_tree_best_per_beam::filter_high_chi_sq_events() {
  for (size_t i = 0; i < event_column_data.size(); ++i) {
    auto pair_key = std::make_pair(event_column_data[i], beam_column_data[i]);
    if (!contains_event_id(pair_key)) {
      update_combo_data(i);
    } else {
      if (event_beam_as_key_map[pair_key].get_chi_sq() >
          chi_sq_column_data[i]) {
        update_combo_data(i);
      }
    }
  }
}

// constructor for compare_hypotheses manager class. initializes two
// hypothesisTrees and the match counter.
compare_hypotheses::compare_hypotheses(
    std::string file_1, std::string tree_1,
    std::vector<Tree_config> alt_hypo_configs, bool match_type)
    : match_by_best_per_beam(match_type),
      matches(0),
      num_hypos(alt_hypo_configs.size()) {
  if (match_by_best_per_beam) {
    tree1 = new hypothesis_tree_best_per_beam(file_1, tree_1,
                                              match_by_best_per_beam);
    alt_hypos.reserve(num_hypos);
    for (Tree_config& config : alt_hypo_configs) {
      alt_hypos.push_back(new hypothesis_tree_best_per_beam(
          config.filename, config.treename, match_by_best_per_beam));
    }
    matched_chi_sqs_by_beam.reserve(alt_hypo_configs.size());
  } else {
    tree1 =
        new hypothesis_tree_best_combo(file_1, tree_1, match_by_best_per_beam);
    alt_hypos.reserve(num_hypos);
    for (Tree_config& config : alt_hypo_configs) {
      alt_hypos.push_back(new hypothesis_tree_best_combo(
          config.filename, config.treename, match_by_best_per_beam));
    }
    matched_chi_sqs.reserve(alt_hypo_configs.size());
  }
}

// load hypothesisTrees' member data from file and cut all high-chisq combos
void compare_hypotheses::prepare_data() {
  tree1->fill_column_vecs();
  tree1->filter_high_chi_sq_events();

  for (hypothesis_tree_base* tree : alt_hypos) {
    tree->fill_column_vecs();
    tree->filter_high_chi_sq_events();
    if (tree->event_column_data.size() == 0) {
      std::cout << "WARNING: Tree " << tree->get_tree_name()
                << " is empty. Did you fill your flat tree?\n";
    }
  }
}

// epsilon is estimated based on smallest difference between different chisq
// values for the same event in dataset
bool compare_hypotheses::chi_sqs_equal(const float& a, const float& b) {
  static constexpr float epsilon = 1e-5f;
  if (std::abs(a - b) < epsilon) {
    return true;
  } else {
    return false;
  }
}

// iterates over tree1's events and checks them against tree2's events.
// events are logged to log_matches.txt and are stored in the matched_chi_sqs
// map.
void compare_hypotheses::find_matches() {
  std::ofstream os;
  if (logging) {
    os.open("log_matches.txt");
  }

  if (!os.good() && logging) {
    std::cerr << "Error: Could not open log file." << std::endl;
    return;
  }

  // match_by_best_per_beam true, match by best combo per beam ID
  if (match_by_best_per_beam) {
    for (hypothesis_tree_base* alt_tree : alt_hypos) {
      std::cout << "Number of unfiltered events in tree1: "
                << tree1->event_column_data.size()
                << " Number of unfiltered events in tree2: "
                << alt_tree->event_column_data.size() << std::endl;
      std::map<std::pair<unsigned long long, unsigned>, float> match_map;
      for (const auto& pair : tree1->event_beam_as_key_map) {
        // get iterator
        auto alt_tree_it = alt_tree->event_beam_as_key_map.find(pair.first);

        // check if key exists
        if (alt_tree_it != alt_tree->event_beam_as_key_map.end()) {
          const combo& alt_combo = alt_tree_it->second;

          // run ID match check
          if (alt_combo.get_run() != pair.second.get_run()) {
            continue;
          }

          // store the match
          match_map[pair.first] = alt_combo.get_chi_sq() / alt_combo.get_ndf();
          matches++;
          if (logging) {
            os << "Event ID: " << (pair.first).first
               << " found in both trees. Run IDs: " << pair.second.get_run()
               << ',' << alt_combo.get_run()
               << " Beam IDs: " << pair.second.get_beam_id() << ','
               << alt_combo.get_beam_id() << '\n';
          }
        }
      }
      // push the map onto compare_hypotheses' vector
      matched_chi_sqs_by_beam.push_back(match_map);
    }
    if (logging) {
      os.close();
    }
    return;
  }

  // match_by_best_per_beam false, match by best overall combo
  for (hypothesis_tree_base* alt_tree : alt_hypos) {
    std::map<unsigned long long, float> match_map;
    for (const auto& pair : tree1->event_as_key_map) {
      // get iterator
      auto alt_tree_it = alt_tree->event_as_key_map.find(pair.first);

      if (alt_tree_it != alt_tree->event_as_key_map.end()) {
        const combo& alt_combo = alt_tree_it->second;

        // run ID match check
        if (alt_combo.get_run() != pair.second.get_run()) {
          continue;
        }

        // store the match
        match_map[pair.first] = alt_combo.get_chi_sq() / alt_combo.get_ndf();
        matches++;

        if (logging) {
          os << "Event ID: " << pair.first
             << " found in both trees. Run IDs: " << pair.second.get_run()
             << ',' << alt_combo.get_run()
             << " Beam IDs: " << pair.second.get_beam_id() << ','
             << alt_combo.get_beam_id() << '\n';
        }
      }
    }
    matched_chi_sqs.push_back(match_map);
  }
  if (logging) {
    os.close();
  }
}

// writes alternative chisq values into new branch. if no match is found,
// placeholder chisq is written instead. if preserve_combos is false (which is
// the default), only the most probable combos from the primary tree and their
// matches are written.
void compare_hypotheses::write_to_file(std::string out_file) {
  if (out_file == "placeholder" || out_file == "") {
    out_file = std::to_string(num_hypos) + "_hypothesesMatched.root";
  }

  float NO_MATCH_INDICATOR = 185100000.0f;

  // write to a computation graph node instead of the actual RDF
  ROOT::RDF::RNode df_node = tree1->df;

  // loop over all alternative hypotheses; add a new branch for each's alt
  // chisqs
  for (size_t i = 0; i < alt_hypos.size(); i++) {
    hypothesis_tree_base* alt_tree = alt_hypos[i];
    auto branchName = alt_tree->get_tree_name() + "_chisq_ndf";

    if (match_by_best_per_beam) {
      auto& matched_chi_sqs_ref = matched_chi_sqs_by_beam[i];
      df_node = df_node.Define(
          branchName,
          [&matched_chi_sqs_ref, NO_MATCH_INDICATOR](unsigned long long event,
                                                     unsigned beam) -> float {
            auto key = std::make_pair(event, beam);
            auto it = matched_chi_sqs_ref.find(key);
            return it != matched_chi_sqs_ref.end() ? it->second
                                                   : NO_MATCH_INDICATOR;
          },
          {"event", "beam_beamid"});
    } else {
      auto& matched_chi_sqs_ref = matched_chi_sqs[i];
      df_node = df_node.Define(branchName,
                               [&matched_chi_sqs_ref, NO_MATCH_INDICATOR](
                                   unsigned long long event) -> float {
                                 auto it = matched_chi_sqs_ref.find(event);
                                 return it != matched_chi_sqs_ref.end()
                                            ? it->second
                                            : NO_MATCH_INDICATOR;
                               },
                               {"event"});
    }
  }

  // preserve only the lowest chisq combo per event ID & beam ID if
  // preserveCombos is false
  if (match_by_best_per_beam) {
    auto& tree1_event_map = tree1->event_beam_as_key_map;
    df_node = df_node.Filter(
        [this, &tree1_event_map](unsigned long long event, unsigned beam,
                                 float kin_chisq) -> bool {
          return preserve_combos ||
                 chi_sqs_equal(kin_chisq,
                               tree1_event_map.at(std::make_pair(event, beam))
                                   .get_chi_sq());
        },
        {"event", "beam_beamid", "kin_chisq"});
  } else {
    // preserve only the lowest chisq combo per event ID if preserveCombos is
    // false
    auto& tree1_event_map = tree1->event_as_key_map;
    df_node = df_node.Filter(
        [this, &tree1_event_map](unsigned long long event,
                                 float kin_chisq) -> bool {
          return preserve_combos ||
                 chi_sqs_equal(kin_chisq,
                               tree1_event_map.at(event).get_chi_sq());
        },
        {"event", "kin_chisq"});
  }

  // process the RNodes and write to file
  df_node.Snapshot("hypothesesMatched", out_file);
}
