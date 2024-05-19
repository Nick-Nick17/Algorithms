#include <algorithm>
#include <bitset>
#include <cstring>
#include <iostream>
#include <stack>
#include <vector>

int LogFloor(int ind) { return __builtin_clzll(1) - __builtin_clzll(ind); }

int LogFloor(size_t ind) { return __builtin_clzll(1) - __builtin_clzll(ind); }

static const size_t kNmbChars = 256;

void PreCalc(int& classes, std::vector<int>& cnt, std::string& str,
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

std::vector<int> GiveSuffixAArray(std::string& str) {
  str += '#';

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

std::vector<int> GiveLCP(const std::string& str,
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

int main() {
  std::string str;
  std::cin >> str;
  std::vector<int> suff_array = GiveSuffixAArray(str);
  std::vector<int> lcp = GiveLCP(str, suff_array);  //  this and next
  std::vector<int> nearest_left = GiveNearestLeft(lcp);
  std::vector<int> nearest_right = GiveNearestRight(lcp);
  SparceTableWithMax max(suff_array);
  SparceTableWithMin min(suff_array);
  int ans = 0;
  std::vector<int> arr(suff_array.size(), 0);
  for (size_t i = 0; i < lcp.size(); ++i) {
    int minimum = min.GiveMinimum(nearest_left[i] + 1, nearest_right[i]);
    int maximum = max.GiveMaximum(nearest_left[i] + 1, nearest_right[i]);
    ans += std::max(0, std::min(lcp[i], maximum - minimum) - arr[i]);
    for (int j = nearest_left[i]; j <= nearest_right[i]; ++j) {
      arr[j] = std::max(arr[j], std::min(lcp[i], maximum - minimum));
    }
  }
  std::cout << ans;
  return 0;
}
