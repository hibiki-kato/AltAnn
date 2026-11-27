#include "../src/altann.hpp"
#include <cassert>
#include <iostream>

int main() {
    const auto fixed = parse_adaptive_k_rule("", 10);
    assert(fixed.budget(1) == 10 && fixed.budget(1000) == 10);
    const auto linear = parse_adaptive_k_rule("linear,4,3", 100);
    const auto exponential = parse_adaptive_k_rule("shifted-exp,5,20,0.1", 40);
    const int exons[] = {1,2,3,5,10,15};
    const int linear_expected[] = {4,7,10,16,31,46};
    const int exponential_expected[] = {5,8,10,15,35,40};
    for (int i = 0; i < 6; ++i) {
        assert(linear.budget(exons[i]) == linear_expected[i]);
        assert(exponential.budget(exons[i]) == exponential_expected[i]);
    }
    assert(exponential.budget(1000000) == 40);
    assert(parse_adaptive_k_rule("shifted-exp,2,0,1e308", 40).budget(1000) == 2);
    assert(parse_adaptive_k_rule("shifted-exp,2,1e308,0", 40).budget(1000) == 2);
    assert(parse_adaptive_k_rule("shifted-exp,2,1,1e308", 40).budget(2) == 40);
    assert(parse_adaptive_k_rule("linear,3,0.5", 40).budget(2) == 4);
    assert(parse_adaptive_k_rule("linear,100,0", 40).budget(1) == 40);
    for (const string &bad : {"linear,0,3", "linear,4,-1", "linear,4,nan", "linear,4,inf",
                              "linear,4,3,", "linear,4,3x", "shifted-exp,5,20,-0.1", "unknown,1,2"}) {
        bool failed = false;
        try { parse_adaptive_k_rule(bad, 40); }
        catch (const runtime_error &) { failed = true; }
        assert(failed);
    }
    std::cout << "Adaptive K examples, clipping and invalid parameters verified\n";
}
