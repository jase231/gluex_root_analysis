#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

struct Tree_config {
  std::string filename;
  std::string treename;
};

struct combo {
 public:
  // getters for member data
  unsigned long long get_event() const { return event; }
  unsigned int get_run() const { return run; }
  unsigned int get_beam_id() const { return beam_beamid; }
  float get_chi_sq() const { return kin_chisq; }
  unsigned get_ndf() const { return kin_ndf; }

  // setters for member data
  void set_chi_sq(float chi_sq) { kin_chisq = chi_sq; }
  void set_event(unsigned long long e) { event = e; }
  void set_run(unsigned int r) { run = r; }
  void set_beam(unsigned int b) { beam_beamid = b; }
  void set_ndf(unsigned n) { kin_ndf = n; }

 private:
  unsigned long long event;
  unsigned int run;
  unsigned int beam_beamid;
  float kin_chisq;
  unsigned kin_ndf;
};

class hypothesis_tree_base {
 public:
  hypothesis_tree_base(std::string file_glob, std::string tree_name,
                       bool match_type);
  virtual ~hypothesis_tree_base() = default;
  // data preperation functions
  virtual void update_combo_data(size_t index) = 0;
  virtual void filter_high_chi_sq_events() = 0;

  bool is_matching_by_beam() const { return match_by_best_per_beam; }
  void set_match_by_beam(bool m) { match_by_best_per_beam = m; }

  bool is_logging() const { return logging; }
  void set_logging(bool l) { logging = l; }

  bool contains_event_id(std::pair<unsigned long long, unsigned>) const;
  void fill_column_vecs();
  std::string get_tree_name() const { return tree_name; }

  // combo maps
  std::map<std::pair<unsigned long long, unsigned>, combo>
      event_beam_as_key_map;
  std::map<unsigned long long, combo> event_as_key_map;

  // RDataFrame
  ROOT::RDataFrame df;

  // column_data vectors
  std::vector<unsigned long long> event_column_data;
  std::vector<unsigned int> run_column_data;
  std::vector<unsigned int> beam_column_data;
  std::vector<float> chi_sq_column_data;
  std::vector<unsigned> ndf_column_data;

 private:
  std::string tree_name;
  bool match_by_best_per_beam;  // whether matching by best combo per beam is
                                // used
  bool logging;
};

class hypothesis_tree_best_combo : public hypothesis_tree_base {
 public:
  hypothesis_tree_best_combo(std::string file_glob, std::string tree_name,
                             bool match_type)
      : hypothesis_tree_base(file_glob, tree_name, match_type) {}
  void update_combo_data(size_t index) override;
  void filter_high_chi_sq_events() override;
};

class hypothesis_tree_best_per_beam : public hypothesis_tree_base {
 public:
  hypothesis_tree_best_per_beam(std::string file_glob, std::string tree_name,
                                bool match_type)
      : hypothesis_tree_base(file_glob, tree_name, match_type) {}
  void update_combo_data(size_t index) override;
  void filter_high_chi_sq_events() override;
};

class compare_hypotheses {
 private:
  hypothesis_tree_base* tree1;
  std::vector<hypothesis_tree_base*> alt_hypos;
  bool logging = false;         // whether logs of matches are written to file
  bool match_by_best_per_beam;  // whether matching by best combo per beam is
                                // used
  bool preserve_combos =
      false;  // whether to keep combos with high chisq in the output file
 public:
  compare_hypotheses(std::string glob1, std::string tree1,
                     std::vector<Tree_config> alt_hypo_configs,
                     bool match_type);

  // calls each tree's data preperation functions
  void prepare_data();

  // performs combo matching between each tree's combo map
  void find_matches();

  // outputs into a new branch in a clone of the primary RDataFrame
  void write_to_file(std::string out_file);

  // float equality function
  bool chi_sqs_equal(const float& a, const float& b);

  // counter for number of matches
  uint matches;
  uint num_hypos;

  // vectors for storing all hypotheses' matches
  std::vector<std::map<std::pair<unsigned long long, unsigned>, float>>
      matched_chi_sqs_by_beam;
  std::vector<std::map<unsigned long long, float>> matched_chi_sqs;

  // helpers and member data setters
  bool is_logging() const { return logging; }
  void set_logging(bool l) {
    logging = l;
    tree1->set_logging(l);
    for (hypothesis_tree_base* tree : alt_hypos) {
      tree->set_logging(l);
    }
  }

  bool is_preserving() const { return preserve_combos; }
  void set_preserving(bool p) { preserve_combos = p; }

  bool is_matching_by_beam() const { return match_by_best_per_beam; }
  void set_match_by_beam(bool m) {
    match_by_best_per_beam = m;
    tree1->set_match_by_beam(m);
    for (hypothesis_tree_base* tree : alt_hypos) {
      tree->set_match_by_beam(m);
    }
  }

  ~compare_hypotheses() {
    delete tree1;
    for (hypothesis_tree_base* tree : alt_hypos) {
      delete tree;
    }
  }
};