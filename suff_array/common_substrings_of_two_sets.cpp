#include <algorithm>
#include <bitset>
#include <cstring>
#include <iostream>
#include <stack>
#include <vector>

static const char kStartChar = 'a';

static int GetIndex(char cc) { return cc - kStartChar + 1; }

int LogFloor(int ind) { return __builtin_clzll(1) - __builtin_clzll(ind); }

int LogFloor(size_t ind) { return __builtin_clzll(1) - __builtin_clzll(ind); }

static const size_t kNmbChars = 2 * 1e6 + 500;

void PreCalc(int& classes, std::vector<int>& cnt, std::vector<int>& str,
             std::vector<int>& cls, std::vector<int>& perm) {
  cnt.assign(std::max(kNmbChars, str.size()), 0);
  for (size_t i = 0; i < str.size(); ++i) {
    ++cnt[static_cast<size_t>(str[i])];
  }
  for (size_t i = 1; i < kNmbChars; ++i) {
    cnt[i] += cnt[i - 1];
  }
  for (int i = static_cast<int>(str.size()) - 1; i >= 0; --i) {
    perm[--cnt[static_cast<size_t>(str[i])]] = i;
  }
  cls[perm[0]] = 0;
  for (size_t i = 1; i < str.size(); ++i) {
    if (str[perm[i]] != str[perm[i - 1]]) {
      ++classes;
    }
    cls[perm[i]] = classes - 1;
  }
}

std::vector<int> GiveSuffixAArray(std::vector<int>& str) {
  str.push_back(0);

  //  Initialization
  int classes = 1;
  std::vector<int> cnt(std::max(kNmbChars, str.size()));
  std::vector<int> perm(str.size());
  std::vector<int> cls(str.size());
  std::vector<int> perm_new(str.size());
  std::vector<int> cls_new(str.size());

  PreCalc(classes, cnt, str, cls, perm);

  //  Sort
  for (size_t height = 0; (1 << height) < str.size(); ++height) {
    for (size_t i = 0; i < str.size(); ++i) {
      perm_new[i] = (static_cast<int>(str.size()) + perm[i] - (1 << height)) %
                    (static_cast<int>(str.size()));
    }
    cnt.assign(classes, 0);
    for (size_t i = 0; i < str.size(); ++i) {
      ++cnt[cls[perm_new[i]]];
    }
    for (int i = 1; i < classes; ++i) {
      cnt[i] += cnt[i - 1];
    }
    for (int i = static_cast<int>(str.size()) - 1; i >= 0; --i) {
      perm[--cnt[cls[perm_new[i]]]] = perm_new[i];
    }
    cls_new[perm[0]] = 0;
    classes = 1;
    for (int i = 1; i < static_cast<int>(str.size()); ++i) {
      bool tmp = cls[(perm[i] + (1 << height)) % str.size()] !=
                 cls[(perm[i - 1] + (1 << height)) % str.size()];
      if (cls[perm[i]] != cls[perm[i - 1]] || tmp) {
        ++classes;
      }
      cls_new[perm[i]] = classes - 1;
    }
    cls.swap(cls_new);
  }

  //  perm.erase(perm.begin());
  //  str.pop_back();
  return perm;
}

std::vector<int> GivePos(const std::vector<int>& suff_array) {
  std::vector<int> pos(suff_array.size());
  for (int i = 0; i < static_cast<int>(suff_array.size()); ++i) {
    pos[suff_array[i]] = i;
  }
  return pos;
}

std::vector<int> GiveLCP(const std::vector<int>& str,
                         const std::vector<int>& suff_array) {
  //  Calculating <pos>
  std::vector<int> pos(suff_array.size());
  std::vector<int> lcp(suff_array.size());
  for (int i = 0; i < static_cast<int>(suff_array.size()); ++i) {
    pos[suff_array[i]] = i;
  }
  int longest = 0;
  int ind = static_cast<int>(suff_array.size());
  for (int word = 0; word < ind; ++word) {
    int index = pos[word];
    if (index + 1 == ind) {
      longest = 0;
      continue;
    }
    int jj = suff_array[index + 1];
    longest = std::max(longest - 1, 0);
    while (word + longest < ind && jj + longest < ind &&
           str[word + longest] == str[jj + longest]) {
      ++longest;
    }
    lcp[index] = longest;
  }
  return lcp;
}

struct SparceTableWithMin {
  std::vector<std::vector<int>> table;
  SparceTableWithMin(const std::vector<int>& lcp) {
    int lg = LogFloor(lcp.size());
    table.resize(lg + 1, std::vector<int>(lcp.size(), -1));
    for (size_t i = 0; i < lcp.size(); ++i) {
      table[0][i] = lcp[i];
    }
    for (int i = 1; i <= lg; ++i) {
      for (int j = 0; j + (1 << i) <= static_cast<int>(lcp.size()); ++j) {
        table[i][j] =
            std::min(table[i - 1][j], table[i - 1][j + (1 << (i - 1))]);
      }
    }
  }

  int GiveMinimum(int left, int right) {
    int index = LogFloor(right - left + 1);
    return std::min(table[index][left], table[index][right - (1 << index) + 1]);
  }
};

struct SparceTableWithMax {
  std::vector<std::vector<int>> table;
  SparceTableWithMax(const std::vector<int>& lcp) {
    int lg = LogFloor(lcp.size());
    table.resize(lg + 1, std::vector<int>(lcp.size(), -1));
    for (size_t i = 0; i < lcp.size(); ++i) {
      table[0][i] = lcp[i];
    }
    for (int i = 1; i <= lg; ++i) {
      for (int j = 0; j + (1 << i) <= static_cast<int>(lcp.size()); ++j) {
        table[i][j] =
            std::max(table[i - 1][j], table[i - 1][j + (1 << (i - 1))]);
      }
    }
  }

  int GiveMaximum(int left, int right) {
    int index = LogFloor(right - left + 1);
    return std::max(table[index][left], table[index][right - (1 << index) + 1]);
  }
};

static const int kMax = 1e9;

std::vector<int> GiveNearestLeft(const std::vector<int>& lcp) {
  std::vector<int> nearest_left(lcp.size(), kMax);
  std::stack<int> stack;
  stack.push(0);
  nearest_left[0] = 0;
  for (size_t i = 1; i < lcp.size(); ++i) {
    while (stack.size() > 1 && lcp[stack.top()] >= lcp[i]) {
      stack.pop();
    }
    nearest_left[i] = stack.top();
    stack.push(i);
  }
  return nearest_left;
}

std::vector<int> GiveNearestRight(const std::vector<int>& lcp) {
  std::vector<int> nearest_right(lcp.size(), kMax);
  std::stack<int> stack;
  stack.push(static_cast<int>(lcp.size()) - 1);
  nearest_right[static_cast<int>(lcp.size()) - 1] =
      static_cast<int>(lcp.size()) - 1;
  for (int i = lcp.size() - 2; i >= 0; --i) {
    while (stack.size() > 1 && lcp[stack.top()] >= lcp[i]) {
      stack.pop();
    }
    nearest_right[i] = stack.top();
    stack.push(i);
  }
  return nearest_right;
}

std::vector<int> GiveSuffNumbers(const std::vector<int>& str,
                                 const std::vector<int>& suff_array,
                                 const std::vector<int>& pos) {
  std::vector<int> suff_numbers(suff_array.size());
  int number = 0;
  for (size_t i = 0; i < str.size(); ++i) {
    if (str[i] > GetIndex('z') || str[i] < GetIndex('a')) {
      suff_numbers[pos[i]] = -1;
      ++number;
      continue;
    }
    suff_numbers[pos[i]] = number;
  }
  return suff_numbers;
}

struct SegTree {
  std::vector<int> tree;
  int size;

  SegTree(std::vector<int>& array) {
    Init(static_cast<int>(array.size()));
    Build(array, 0, 0, size);
  }

  void Init(int vert_number) {
    size = 1;
    while (size < vert_number) {
      size *= 2;
    }
    tree.assign(2 * size - 1, 0);
  }

  void Build(std::vector<int>& array, int x_now, int lx, int rx) {
    if (rx - lx == 1) {
      if (lx < static_cast<int>(array.size())) {
        tree[x_now] = array[lx];
      }
    } else {
      int mid = (lx + rx) / 2;
      Build(array, 2 * x_now + 1, lx, mid);
      Build(array, 2 * x_now + 2, mid, rx);
      tree[x_now] = tree[2 * x_now + 1] + tree[2 * x_now + 2];
    }
  }

  void Set(int index, int vertex, int x_now, int lx, int rx) {
    if (rx - lx == 1) {
      tree[x_now] = vertex;
      return;
    }
    int mid = (lx + rx) / 2;
    if (index < mid) {
      Set(index, vertex, 2 * x_now + 1, lx, mid);
    } else {
      Set(index, vertex, 2 * x_now + 2, mid, rx);
    }
    tree[x_now] = tree[2 * x_now + 1] + tree[2 * x_now + 2];
  }

  void Set(int index, int vertex) { Set(index, vertex, 0, 0, size); }

  int Sum(int left, int right, int x_now, int lx, int rx) {
    if (left >= rx || lx >= right) {
      return 0;
    }
    if (lx >= left && rx <= right) {
      return tree[x_now];
    }
    int mid = (lx + rx) / 2;
    return Sum(left, right, 2 * x_now + 1, lx, mid) +
           Sum(left, right, 2 * x_now + 2, mid, rx);
  }

  int Sum(int left, int right) { return Sum(left, right, 0, 0, size); }
};

static const int kMaxL = 2 * 1e6;

std::vector<int> GiveNumberOfDiff(std::vector<int>& array,
                                  std::vector<std::pair<int, int>>& segments) {
  std::vector<int> ans(segments.size(), -1);
  std::vector<std::pair<int, int>> r_seg;
  r_seg.reserve(segments.size());
  for (int i = 0; i < static_cast<int>(segments.size()); ++i) {
    r_seg.push_back({segments[i].second, i});
  }
  std::sort(r_seg.begin(), r_seg.end());
  std::vector<int> sta(array.size(), 0);
  SegTree states(sta);
  std::vector<int> pointer_before(kMaxL, -1);
  std::vector<int> next(array.size(), -1);
  //  массив, в котором храню куда число
  int segment_now = 0;
  for (int i = 0; i < static_cast<int>(array.size()); ++i) {
    if (array[i] != -1) {
      int latest = pointer_before[array[i]];
      if (latest != -1) {
        states.Set(latest, 0);
      }
      states.Set(i, 1);
      pointer_before[array[i]] = i;
    }
    while (segment_now < static_cast<int>(segments.size()) &&
           r_seg[segment_now].first == i) {
      int seg_to_find = r_seg[segment_now].second;
      ans[seg_to_find] = states.Sum(0, i + 1);
      ans[seg_to_find] -= states.Sum(0, segments[seg_to_find].first);
      ++segment_now;
    }
  }
  return ans;
}

std::vector<std::pair<int, int>> GiveSegments(std::vector<int>& nearest_left,
                                              std::vector<int>& nearest_right) {
  std::vector<std::pair<int, int>> res;
  res.reserve(nearest_left.size());
  for (size_t i = 0; i < nearest_left.size(); ++i) {
    res.push_back({nearest_left[i] + 1, nearest_right[i]});
  }
  return res;
}

void FindAns(std::vector<int> str, int ll) {
  std::vector<int> suff_array = GiveSuffixAArray(str);
  std::vector<int> lcp = GiveLCP(str, suff_array);
  std::vector<int> pos = GivePos(suff_array);
  std::vector<int> suff_numbers = GiveSuffNumbers(str, suff_array, pos);
  SparceTableWithMin min(lcp);
  std::vector<int> nearest_left = GiveNearestLeft(lcp);
  std::vector<int> nearest_right = GiveNearestRight(lcp);
  std::vector<std::pair<int, int>> segments =
      GiveSegments(nearest_left, nearest_right);
  std::vector<int> diff_numbers = GiveNumberOfDiff(suff_numbers, segments);
  std::vector<int> ans(ll, 0);
  for (size_t i = 0; i < suff_array.size(); ++i) {
    int index = diff_numbers[i];
    ans[index - 1] = std::max(ans[index - 1], lcp[i]);
  }
  for (int i = ll - 1; i >= 1; --i) {
    ans[i - 1] = std::max(ans[i - 1], ans[i]);
  }
  for (int i = 1; i < ll; ++i) {
    std::cout << ans[i] << '\n';
  }
}

int main() {
  int str_number;
  std::cin >> str_number;
  std::vector<int> string;

  for (int i = 0; i < str_number; ++i) {
    if (i != 0) {
      string.push_back(GetIndex('z') + 1 + i);
    }
    std::string tmp;
    std::cin >> tmp;
    for (size_t j = 0; j < tmp.size(); ++j) {
      string.push_back(GetIndex(tmp[j]));
    }
  }
  FindAns(string, str_number);
}
