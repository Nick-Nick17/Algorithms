#include <algorithm>
#include <iostream>
#include <vector>

bool Bit(int mask, int pos) { return ((mask >> pos) & 1) != 0; }

bool Prove(int m, int mask, std::vector<std::vector<int>>& left_change,
           std::vector<std::vector<int>>& right_change) {
  for (auto number : left_change[m]) {
    if (Bit(mask, number)) {
      return false;
    }
  }
  for (auto number : right_change[m]) {
    if (!Bit(mask, number)) {
      return false;
    }
  }
  return true;
}

void Solution() {
  const long long kMod = 1e9 + 7;
  int n, m;
  std::cin >> n >> m;
  int full = 0;
  int tried_left = 0, tried_right = 0;
  for (int i = 0; i < n; ++i) {
    full += (1 << i);
    if (i % 2 == 0) {
      tried_left += (1 << i);
    } else {
      tried_right += (1 << i);
    }
  }

  std::vector<std::vector<int>> left_change(m);
  std::vector<std::vector<int>> right_change(m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      char c;
      std::cin >> c;
      if (c == '+') {
        left_change[j].push_back(i);
      }
      if (c == '-') {
        right_change[j].push_back(i);
      }
    }
  }

  std::vector<std::vector<long long>> dp((1 << n),
                                         std::vector<long long>(m + 1));
  for (int i = 0; i < m; ++i) {
    for (int mask = 0; mask < (1 << n); ++mask) {
      dp[mask][i] = 0;
      if (Prove(i, mask, left_change, right_change)) {
        if (i == 0) {
          dp[mask][i] = 1;
        } else {
          int left = (~mask) & full;
          dp[mask][i] = dp[left][i - 1];
          if (mask == tried_right || mask == tried_left) {
            dp[mask][i] += dp[mask][i - 1];
          }
          dp[mask][i] %= kMod;
        }
      }
    }
  }
  long long ans = 0;
  for (int i = 0; i < (1 << n); ++i) {
    ans = (ans + dp[i][m - 1]) % kMod;
  }
  std::cout << ans;
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(0);
  std::cout.tie(0);
  Solution();
  return 0;
}
