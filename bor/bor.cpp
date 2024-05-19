#include <algorithm>
#include <bitset>
#include <cstring>
#include <iostream>
#include <queue>
#include <vector>

static const size_t kLibrary = 2;

struct Node {
  bool is_terminate;
  int link;
  int to[kLibrary];

  int compresed = -1;
  int dp = 0;
  std::vector<int> terminate_number;

  Node() : is_terminate(false), link(-1) { std::memset(to, -1, sizeof(to)); }
};

struct Bor {
  int last;
  std::string text;
  std::vector<Node> data;
  std::vector<int> answer;
  std::vector<std::vector<int>> all_in;

  std::vector<std::string> strings;

  Bor() {
    last = 1;
    size_t number;
    std::cin >> number;
    strings.resize(number);
    answer.resize(number, 0);
    all_in.resize(number);
    data.push_back(Node());
    for (size_t i = 0; i < number; ++i) {
      std::cin >> strings[i];
      AddString(i);
    }
  }

  static size_t GetIndex(char in) { return in - '0'; }

  void AddString(size_t str) { AddString(str, 0, 0); }

  void AddString(size_t str, size_t index, int node_number) {
    if (index == strings[str].size()) {
      data[node_number].is_terminate = true;
      data[node_number].terminate_number.push_back(str);
      ++data[node_number].dp;
      return;
    }
    size_t tmp = GetIndex(strings[str][index]);
    if (data[node_number].to[tmp] == -1) {
      data.push_back(Node());
      data[node_number].to[tmp] = last;
      ++last;
    }
    AddString(str, index + 1, data[node_number].to[tmp]);
  }

  bool DFS(int vertex, std::vector<int>& used) {
    used[vertex] = 0;
    for (size_t i = 0; i < kLibrary; ++i) {
      int next = data[vertex].to[i];
      if (data[next].is_terminate || data[next].dp != 0) {
        continue;
      }
      if (used[next] == 1) {
        continue;
      }
      if (used[next] == 0 || DFS(next, used)) {
        return true;
      }
    }
    used[vertex] = 1;
    return false;
  }

  void FindAns() {
    std::vector<int> used(data.size(), -1);
    if (data[0].to[0] == 0 || data[0].to[1] == 0) {
      std::cout << "ТАК";
      return;
    }
    if (DFS(0, used)) {
      std::cout << "TAK";
    } else {
      std::cout << "NIE";
    }
  }
};

void Inside(Bor& bor, int next, int vertex) {
  if (bor.data[bor.data[next].link].is_terminate) {
    bor.data[next].compresed = bor.data[next].link;
  } else {
    bor.data[next].compresed = bor.data[bor.data[next].link].compresed;
  }
  bor.data[next].dp += bor.data[bor.data[next].link].dp;
  bor.data[next].dp += bor.data[vertex].dp;
}

void AhoCorasik(Bor& bor) {
  bor.data[0].link = 0;
  for (size_t cc = 0; cc < kLibrary; ++cc) {
    if (bor.data[0].to[cc] != -1) {
      continue;
    }
    bor.data[0].to[cc] = 0;
  }
  std::queue<int> qq;
  qq.push(0);
  while (!qq.empty()) {
    int vv = qq.front();
    qq.pop();
    for (size_t i = 0; i < kLibrary; ++i) {
      int next = bor.data[vv].to[i];
      if (bor.data[next].link != -1) {
        continue;
      }
      bor.data[next].link = (vv == 0 ? 0 : bor.data[bor.data[vv].link].to[i]);
      for (size_t cc = 0; cc < kLibrary; ++cc) {
        if (bor.data[next].to[cc] != -1) {
          continue;
        }
        bor.data[next].to[cc] = bor.data[bor.data[next].link].to[cc];
      }
      Inside(bor, next, vv);
      qq.push(next);
    }
  }
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(0);
  std::cout.tie(0);

  Bor bor;
  AhoCorasik(bor);
  bor.FindAns();
  return 0;
}
