#ifndef TWOSTATESEQ_HPP
#define TWISTATESEQ_HPP

#include <string>
#include <vector>
#include <random>

class TwoStateSeq {
public:
  void set_params(size_t n, double p0, double p1);
  void simulate();
  void count_s3(std::vector<size_t> &freq) const; 
  void mutation(const size_t pos, cosnt size_t s3_pattern,
                std::vector<size_t> &freq);
  
private:
  std::vector<double> params;  // two self-transition probabilities
  size_t n_site; 
  std::vector<std::tuple<size_t, size_t> > fg_intervals;
};


TwoStateSeq::set_params(size_t n, double p0, double p1) {
  n_site = n;
  assert (p0 > 0 && p1 > 0 && p0 < 1 && p1 < 1);
  params.resize(2);
  params[0] = p0;
  params[1] = p1;
}

TwoStateSeq::simulate() {
  // Assume first and last sites are fixed in the background state  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::geometric_distribution<size_t> fg_geo(1.0 - params[0]); // non-negative integer
  std::geometric_distribution<size_t> bg_geo(1.0 - params[1]);  
  size_t start = bg_geo(gen) + 1;
  while (start < n_site - 1) {
    size_t end = std::min(start + fg_geo(gen), n_site - 2);
    fg_intervals.push_back(std::make_tuple(start, end));
    start = end + 1 + bg_geo(gen);
  }
}


TwoStateSeq::count_s3(std::vector<size_t> &freq) {
  freq.resize(8, 0);
  // 111
  freq[7] += std::max(0, fg_intervals[0].first - 2); 
  freq[7] += std::max(0, n_site - 3 - fg_intervals.back().second);
  for (size_t i = 0; i < fg_intervals.size(); ++i) {
    size_t l = fg_intervals[i].first;
    size_t r = fg_intervals[i].second;
    freq[0] += std::max(0, r - l - 1);     // 000
    if (i > 0)
      freq[7] += std::max(0, l - fg_intervals[i-1].second -3);     // 111
    freq[2] += (size_t)(i > 0 && l == fg_intervals[i-1].second + 2);     // 010
    freq[5] += (size_t)(l == r); // 101
    freq[1] += (size_t)(l < r); // 001
    freq[4] += (size_t)(l < r); // 100
    // 011
    freq[3] += (size_t)(i < fg_intervals.size() - 1 && r + 2 < fg_intervals[i+1].first);
    freq[3] += (size_t)(i == fg_intervals.size() - 1 && r + 2 < n_site);
    //110
    freq[6] += (size_t)((i == 0 && l > 1)|(i > 0 && l > fg_intervals[i-1].second + 2));
  }
}


