#include <algorithm>
#include <bitset>
#include <cstring>
#include <iostream>
#include <list>
#include <vector>

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
  return perm;
}

std::string MTF(std::string& str) {
  std::string res;
  res.reserve(str.size());
  std::list<char> alphabet_for_mtf;
  for (int i = 0; i <= 'z' - 'a'; ++i) {
    alphabet_for_mtf.push_back(static_cast<char>('a' + i));
  }
  for (size_t i = 0; i < str.size(); ++i) {
    auto it = alphabet_for_mtf.begin();
    int index = 0;
    while (*it != str[i]) {
      ++index;
      ++it;
    }
    res.push_back(index + 'a');
    alphabet_for_mtf.erase(it);
    alphabet_for_mtf.push_front(str[i]);
  }
  return res;
}

std::string DeMTF(std::string& str) {
  std::string res;
  res.reserve(str.size());
  std::list<char> alphabet_for_mtf;
  for (int i = 0; i <= 'z' - 'a'; ++i) {
    alphabet_for_mtf.push_back(static_cast<char>('a' + i));
  }
  for (size_t i = 0; i < str.size(); ++i) {
    auto it = alphabet_for_mtf.begin();
    int index = str[i] - 'a';
    while (index > 0) {
      --index;
      ++it;
    }
    res.push_back(*it);
    alphabet_for_mtf.erase(it);
    alphabet_for_mtf.push_front(res[res.size() - 1]);
  }
  return res;
}

void RLE(std::string& str) {
  char symb = str[0];
  int index = 1;
  for (size_t i = 1; i < str.size(); ++i) {
    if (symb == str[i]) {
      ++index;
    } else {
      std::cout << symb << index;
      index = 1;
      symb = str[i];
    }
  }
  std::cout << symb << index;
}

void Code() {
  std::string str;
  std::cin >> str;
  int numb = static_cast<int>(str.size());
  std::vector<int> my_suffix_array = GiveSuffixAArray(str);
  std::string bwt_res;
  bwt_res.reserve(str.size());
  for (size_t i = 0; i < str.size(); ++i) {
    bwt_res.push_back(str[(numb + my_suffix_array[i] - 1) % numb]);
  }
  std::string mtf_res = MTF(bwt_res);
  RLE(mtf_res);
}

std::string DecodingRLE(std::string& str) {
  std::string res;
  for (size_t i = 0; i < str.size(); ++i) {
    char symb = str[i];
    int number = 0;
    while (i + 1 < str.size() && '0' <= str[i + 1] && str[i + 1] <= '9') {
      number = number * (3 + 3 + 3 + 1) + str[i + 1] - '0';
      ++i;
    }
    for (int j = 0; j < number; ++j) {
      res += symb;
    }
  }
  return res;
}

void DecodingBWT(std::string& str) {
  int perm_nmb;
  std::cin >> perm_nmb;
  std::string old_str = str;
  std::sort(str.begin(), str.end());
  std::vector<int> res(str.size());
  std::vector<std::vector<int>> arr_lists('z' - 'a' + 1);
  for (int i = 0; i < static_cast<int>(str.size()); ++i) {
    arr_lists[old_str[i] - 'a'].push_back(i);
  }
  std::vector<size_t> sdvig('z' - 'a' + 1, 0);
  for (size_t i = 0; i < str.size(); ++i) {
    res[i] = arr_lists[str[i] - 'a'][sdvig[str[i] - 'a']];
    ++sdvig[str[i] - 'a'];
  }
  std::vector<char> answer(str.size());
  for (size_t i = 0; i < str.size(); ++i) {
    perm_nmb = res[perm_nmb];
    answer[str.size() - i - 1] = old_str[perm_nmb];
  }
  for (int i = static_cast<int>(str.size()) - 1; i >= 0; --i) {
    std::cout << answer[i];
  }
}

void DeCode() {
  std::string str;
  std::cin >> str;
  std::string res = DecodingRLE(str);
  std::string from_mtf = DeMTF(res);
  DecodingBWT(from_mtf);
}

int main() {
  int type;
  std::cin >> type;
  if (type == 1) {
    Code();
  } else {
    DeCode();
  }
  return 0;
}
