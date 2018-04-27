#ifndef TRIPLETPATTERN_HPP
#define TRIPLETPATTERN_HPP

#include <vector>
#include <random>
#include <algorithm>

class PatSeq {
public:
  PatSeq(const std::vector<bool> &seq);
  void get_all_context_freq(std::vector<size_t> &freq) const;
  size_t get_context(const size_t pos) const;
  size_t get_context_freq(const size_t context) const;
  size_t random_mutate(const size_t context, std::mt19937 &gen);
  void mutate(const size_t pos);
  void to_seq(std::vector<bool> &seq) const;
  void to_domain_sizes(std::vector<size_t> &domain_sizes) const;
  size_t get_size() {return pos_by_pat.size();};

private:
  /* positions in binary-state sequence
     organized by triplet pattern */
  std::vector<size_t> pos_by_pat;
  /* indices in pos_by_pat for
     positions in binary-state sequence*/
  std::vector<size_t> idx_in_pat;
  std::vector<size_t> cum_pat_freq;
  bool start_state;
  bool end_state;
  void single_update(const size_t pos, const size_t context,
                     const size_t to_context);
};


#endif
