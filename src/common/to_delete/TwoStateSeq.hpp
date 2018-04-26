#ifndef TWOSTATESEQ_HPP
#define TWOSTATESEQ_HPP

#include <algorithm>    // std::max

#include <string>
#include <vector>
#include <random>

class TwoStateSeq {
public:
  void set_params(size_t n, double p0, double p1);
  void simulate();
  void count_s3(std::vector<size_t> &freq) const;
  void mutate(const size_t pos, size_t &context);
  TwoStateSeq(const std::vector<bool> &seq);

private:
  std::vector<double> params;  // two self-transition probabilities
  size_t n_site;
  std::vector<std::pair<size_t, size_t> > fg_intervals; // both ends inclusive
};

void TwoStateSeq::set_params(size_t n, double p0, double p1) {
  n_site = n;
  assert (p0 > 0 && p1 > 0 && p0 < 1 && p1 < 1);
  params.resize(2);
  params[0] = p0;
  params[1] = p1;
}

void TwoStateSeq::simulate() {
  // Assume first and last sites are fixed in the background state
  std::random_device rd;
  std::mt19937 gen(rd());
  std::geometric_distribution<size_t> fg_geo(1.0 - params[0]); // non-negative integer
  std::geometric_distribution<size_t> bg_geo(1.0 - params[1]);
  size_t start = bg_geo(gen) + 1;
  while (start < n_site - 1) {
    size_t end = std::min(start + fg_geo(gen), n_site - 2);
    fg_intervals.push_back(std::make_pair(start, end));
    start = end + 1 + bg_geo(gen);
  }
}

void TwoStateSeq::count_s3(std::vector<size_t> &freq) const {
  freq.resize(8, 0);
  if (fg_intervals.size() == 0) {
    freq[7] = n_site - 2; 
  } else {
    size_t l = fg_intervals[0].first;
    size_t r = fg_intervals[0].second;
    size_t prev_r = 0;
    // 111
    if (l > 2) freq[7] += l - 2; // 111 before first interval
    if (l > 0 && l < r) ++freq[4];  // 100
    if (l > 1) ++freq[6]; // 110 (0 in first interval)
    for (size_t i = 0; i < fg_intervals.size(); ++i) {
      l = fg_intervals[i].first;
      r = fg_intervals[i].second;
      if (r > l + 1)
        freq[0] += r - l - 1;     // 000
      if (i > 0) {
        if (l > prev_r + 3 )
          freq[7] += l - prev_r -3; // 111 ( before this interval)
        if (l == prev_r + 2)
          ++ freq[2]; // 010
        if (l > prev_r +2) {
          ++ freq[3]; // 011 (0 in previous interval)
          ++ freq[6]; // 110 (0 in current interval)
        }        
      } 
      if (l == r) ++freq[5]; // 101
      if (l < r) {
        if (r < n_site - 1) ++freq[1];  // 001
        ++freq[4];  // 100
      }
      prev_r = r;
    }
    if (r+3 < n_site) freq[7] += n_site - r - 3; // after last interval
  }
}


void TwoStateSeq::mutate(const size_t pos, size_t &context) {
  assert(pos > 0 and pos < n_site - 1);
  // Locate pos relative to fg_intervals
  size_t i = 0;
  while (i < fg_intervals.size() and fg_intervals[i].second < pos) ++i;
  if (i < fg_intervals.size()) {
    if (pos < fg_intervals[i].first -1) {
      // 111; insert
      context = 7;
      fg_intervals.insert(fg_intervals.begin()+i, std::make_pair(pos, pos));
    } else if (pos == fg_intervals[i].first && pos == fg_intervals[i].second) {
      // 101; delete
      context = 5;
      fg_intervals.erase(fg_intervals.begin()+i);
    } else if (pos == fg_intervals[i].first && pos < fg_intervals[i].second) {
      // 100; contraction
      context = 4;
      ++ fg_intervals[i].first;
    } else if (pos == fg_intervals[i].second &&
               fg_intervals[i].second > fg_intervals[i].first) {
      // 001; contraction
      context = 1;
      -- fg_intervals[i].second;
    } else if (pos > fg_intervals[i].first &&
               pos < fg_intervals[i].second) {
      // 000; split
      context = 0;
      fg_intervals.insert(fg_intervals.begin()+i, fg_intervals[i]);
      fg_intervals[i].second = pos - 1;
      fg_intervals[i+1].first = pos + 1;
    } else if (i > 0 && pos == fg_intervals[i-1].second + 1 &&
               pos == fg_intervals[i].first -1) {
      // 010;  merge
      context = 2;
      fg_intervals[i].first = fg_intervals[i-1].first;
      fg_intervals.erase(fg_intervals.begin()+i-1);
    } else if (i > 0 && pos == fg_intervals[i-1].second + 1 &&
               pos < fg_intervals[i].first -1) {
      // 011; extension
      context = 3;
      ++ fg_intervals[i-1].second;
    } else {  // the last possibility: i == 0 and pos == fg_intervals[i].first -1
      // 110; extension
      context = 6;
      --fg_intervals[i].first;
    }
  } else {
    if (pos == fg_intervals.back().second + 1) {
      // 011; extension
      context = 3;
      ++ fg_intervals.back().second;
    } else {
      // 111; insert
      context = 7;
      fg_intervals.push_back(std::make_pair(pos, pos));
    }
  }
}


TwoStateSeq::TwoStateSeq(const std::vector<bool> &seq) {
  // params won't be initialized
  fg_intervals.clear();
  n_site = seq.size();
  bool in_domain = false;
  size_t start = 0;
  size_t end = 0;
  for (size_t i = 0; i < seq.size(); ++i) {
    if (!seq[i] && !in_domain) {
      start = i;
      in_domain = true;
    } else if (seq[i] && in_domain) {
      end = i-1;
      fg_intervals.push_back(std::make_pair(start, end));
      in_domain = false;
    }
  }

  if (in_domain) {
    fg_intervals.push_back(std::make_pair(start, n_site-1));
  }
}

#endif
