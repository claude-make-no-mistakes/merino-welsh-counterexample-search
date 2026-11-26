#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Count = std::int64_t;
using Wide = unsigned __int128;
constexpr int edge_limit = 20;

struct Bundle {
    int u, v, copies;
};

struct Values {
    Count acyclic;
    Count cyclic;
    Count trees;
};

class Decorations {
public:
    Decorations(int n, const std::vector<Bundle>& edges)
        : n(n), edges(edges), a(1U << edges.size()), c(a.size()), t(a.size()) {
        std::iota(parent.begin(), parent.end(), 0);
        component_size.fill(1);
        ranks(0, 0, 0, n);
        transform();
    }

    std::size_t size() const { return a.size(); }
    Values at(unsigned mask) const { return {a[mask], c[mask], t[mask]}; }

private:
    int n;
    const std::vector<Bundle>& edges;
    std::array<int, edge_limit> parent, component_size;
    std::vector<Count> a, c, t;

    int root(int vertex) const {
        while (parent[vertex] != vertex) {
            vertex = parent[vertex];
        }
        return vertex;
    }

    void ranks(unsigned bit, unsigned subset, int selected, int components) {
        if (bit == edges.size()) {
            a[subset] = ((n + components) % 2 == 0) ? 1 : -1;
            c[subset] = (components % 2 == 1) ? 1 : -1;
            t[subset] = (components == 1 && selected == n - 1);
            return;
        }
        ranks(bit + 1, subset, selected, components);
        const auto [u, v, copies] = edges[bit];
        int left = root(u), right = root(v);
        const unsigned include = subset | (1U << bit);
        if (left == right) {
            ranks(bit + 1, include, selected + 1, components);
            return;
        }
        if (component_size[left] < component_size[right]) {
            std::swap(left, right);
        }
        parent[right] = left;
        component_size[left] += component_size[right];
        ranks(bit + 1, include, selected + 1, components - 1);
        component_size[left] -= component_size[right];
        parent[right] = right;
    }

    void transform() {
        for (unsigned bit = 0; bit < edges.size(); ++bit) {
            const unsigned stride = 1U << bit;
            const int copies = edges[bit].copies;
            const Count power2 = Count{1} << copies;
            Count power3 = 1;
            for (int i = 1; i < copies; ++i) {
                power3 *= 3;
            }
            // Row 0 has one doubled edge. Row 1 subdivides the entire bundle.
            for (unsigned block = 0; block < a.size(); block += 2 * stride) {
                for (unsigned lo = block; lo < block + stride; ++lo) {
                    const unsigned hi = lo + stride;
                    const Count a0 = a[lo], a1 = a[hi];
                    const Count c0 = c[lo], c1 = c[hi];
                    const Count t0 = t[lo], t1 = t[hi];
                    a[lo] = power3 * (a0 - a1);
                    a[hi] = 3 * power3 * a0 + (power2 - 3 * power3) * a1;
                    c[lo] = c0 + (2 * power2 - 1) * c1;
                    c[hi] = c0 + (power2 - 1) * c1;
                    t[lo] = (power2 / 2) * t0 + ((copies + 3) * power2 / 4) * t1;
                    t[hi] = power2 * t0 + (copies * (power2 / 2)) * t1;
                }
            }
        }
    }
};

bool less_fraction(Wide a, Wide b, Wide c, Wide d) {
    // Comparing continued fractions avoids products of four counts.
    bool flip = false;
    while (true) {
        if (a / b != c / d) {
            return flip ? a / b > c / d : a / b < c / d;
        }
        a %= b;
        c %= d;
        if (a == 0 || c == 0) {
            return flip ? a > c : a < c;
        }
        std::swap(a, b);
        std::swap(c, d);
        flip = !flip;
    }
}

bool series_parallel(int n, const std::vector<Bundle>& edges) {
    // Simplify parallel edges, then complete and remove degree-two vertices.
    std::array<unsigned, edge_limit> graph{};
    for (const auto& [u, v, copies] : edges) {
        graph[u] |= 1U << v;
        graph[v] |= 1U << u;
    }
    unsigned active = (1U << n) - 1;
    while (active) {
        int u = 0;
        while (u < n && (((active >> u) & 1U) == 0 || __builtin_popcount(graph[u]) > 2)) {
            ++u;
        }
        if (u == n) {
            return false;
        }
        unsigned neighbors = graph[u];
        while (neighbors) {
            const int v = __builtin_ctz(neighbors);
            neighbors &= neighbors - 1;
            graph[v] = (graph[v] & ~(1U << u)) | (graph[u] & ~(1U << v));
        }
        active &= ~(1U << u);
    }
    return true;
}

struct Minimum {
    Values counts{};
    unsigned mask = 0;
    std::vector<Bundle> edges;
    Wide numerator = 0;
    Wide denominator = 0;

    void consider(const Values& values, unsigned decoration, const std::vector<Bundle>& core) {
        const Wide a = static_cast<Wide>(values.acyclic) * values.cyclic;
        const Wide b = static_cast<Wide>(values.trees) * values.trees;
        // All four-count products fit in 128 bits through m=15.
        if (denominator != 0 &&
            !(a < (Wide{1} << 64) && b < (Wide{1} << 64) &&
              numerator < (Wide{1} << 64) && denominator < (Wide{1} << 64)
                ? a * denominator < numerator * b
                : less_fraction(a, b, numerator, denominator))) {
            return;
        }
        counts = values;
        mask = decoration;
        edges = core;
        numerator = a;
        denominator = b;
    }

    void print() const {
        if (denominator == 0) {
            std::cout << "null";
            return;
        }
        std::cout << "{\"series_mask\":" << mask << ",\"acyclic\":" << counts.acyclic
                  << ",\"cyclic\":" << counts.cyclic << ",\"trees\":" << counts.trees
                  << ",\"bundles\":[";
        for (unsigned i = 0; i < edges.size(); ++i) {
            if (i != 0) {
                std::cout << ',';
            }
            std::cout << '[' << edges[i].u << ',' << edges[i].v << ',' << edges[i].copies << ']';
        }
        std::cout << "]}";
    }
};

struct Result {
    std::uint64_t cores = 0, series_parallel_cores = 0, cases = 0, equality = 0, failures = 0;
    Minimum minimum, counterexample;
};

void scan(bool full_table) {
    std::array<std::array<Result, edge_limit + 1>, edge_limit + 1> results;
    int n, unique_edges;
    while (std::cin >> n) {
        if (!(std::cin >> unique_edges) || n < 2 || n > edge_limit
            || unique_edges < 1 || unique_edges > edge_limit) {
            throw std::runtime_error("invalid graph header");
        }
        std::vector<Bundle> edges;
        unsigned m = 0;
        for (int i = 0; i < unique_edges; ++i) {
            int u, v, copies;
            if (!(std::cin >> u >> v >> copies) || u < 0 || u >= v || v >= n || copies < 1
                || copies > edge_limit - static_cast<int>(m)) {
                throw std::runtime_error("invalid edge row or too many edges");
            }
            edges.push_back({u, v, copies});
            m += copies;
        }
        auto& result = results[n][m];
        ++result.cores;
        const bool sp = series_parallel(n, edges);
        result.series_parallel_cores += sp;
        if (sp && !full_table) {
            continue;
        }
        const Decorations table(n, edges);
        result.cases += table.size();
        for (unsigned mask = 0; mask < table.size(); ++mask) {
            const auto values = table.at(mask);
            if (values.acyclic <= 0 || values.cyclic <= 0 || values.trees <= 0) {
                throw std::runtime_error("nonpositive evaluation for a loopless bridgeless core");
            }
            const Wide product = static_cast<Wide>(values.acyclic) * values.cyclic;
            const Wide square = static_cast<Wide>(values.trees) * values.trees;
            result.equality += product == square;
            result.failures += product < square;
            result.minimum.consider(values, mask, edges);
            if (product < square) {
                result.counterexample.consider(values, mask, edges);
            }
            if (full_table) {
                std::cout << '[' << mask << ',' << values.acyclic << ','
                          << values.cyclic << ',' << values.trees << "]\n";
            }
        }
    }
    if (!std::cin.eof()) {
        throw std::runtime_error("invalid graph stream");
    }
    if (full_table) {
        return;
    }
    for (int v = 2; v <= edge_limit; ++v) {
        for (int m = 2; m <= edge_limit; ++m) {
            const auto& row = results[v][m];
            if (row.cores == 0) {
                continue;
            }
            std::cout << "{\"vertices\":" << v << ",\"core_edges\":" << m
                      << ",\"cores\":" << row.cores << ",\"series_parallel_cores\":" << row.series_parallel_cores
                      << ",\"decorations\":" << row.cases << ",\"equalities\":"
                      << row.equality << ",\"failures\":" << row.failures << ",\"minimum\":";
            row.minimum.print();
            std::cout << ",\"counterexample\":";
            row.counterexample.print();
            std::cout << "}\n";
        }
    }
}

int main(int argc, char** argv) {
    const bool full_table = argc == 2 && std::string(argv[1]) == "--table";
    if (argc > 1 && !full_table) {
        std::cerr << "usage: search [--table]\n";
        return 2;
    }
    std::ios::sync_with_stdio(false);
    try {
        scan(full_table);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
