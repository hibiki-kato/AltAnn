#include "altann.hpp"

namespace {
double nonnegative_number(const string &text) {
    size_t used = 0;
    double value;
    try { value = stod(text, &used); }
    catch (const exception &) { throw runtime_error("Invalid adaptive K number: " + text); }
    if (used != text.size() || !isfinite(value) || value < 0)
        throw runtime_error("Adaptive K parameters must be finite and nonnegative");
    return value;
}
}

AdaptiveKRule parse_adaptive_k_rule(const string &spec, int cap) {
    if (cap < 1 || cap > 32767) throw runtime_error("Adaptive K cap must be between 1 and 32767");
    AdaptiveKRule rule;
    rule.cap = cap;
    rule.base = cap;
    if (spec.empty()) return rule;
    vector<string> fields;
    istringstream input(spec);
    string field;
    while (getline(input, field, ',')) fields.push_back(field);
    if (spec.back() == ',' || fields.empty()) throw runtime_error("Invalid adaptive K specification");
    if (fields[0] == "linear" && fields.size() == 3)
        rule.family = AdaptiveKRule::Family::Linear;
    else if (fields[0] == "shifted-exp" && fields.size() == 4)
        rule.family = AdaptiveKRule::Family::ShiftedExponential;
    else throw runtime_error("Adaptive K must be linear,BASE,SLOPE or shifted-exp,BASE,SCALE,RATE");
    rule.base = nonnegative_number(fields[1]);
    if (rule.base < 1) throw runtime_error("Adaptive K base must be at least 1");
    rule.scale = nonnegative_number(fields[2]);
    if (rule.family == AdaptiveKRule::Family::ShiftedExponential)
        rule.rate = nonnegative_number(fields[3]);
    return rule;
}

int AdaptiveKRule::budget(size_t exon_count) const {
    if (exon_count < 1) throw runtime_error("Adaptive K requires a nonempty CDS chain");
    if (family == Family::Fixed) return cap;
    if (base >= cap) return cap;
    const long double steps = exon_count - 1;
    long double value = base;
    if (family == Family::Linear) {
        value += static_cast<long double>(scale) * steps;
    } else if (scale > 0 && rate > 0 && steps > 0) {
        // Test the cap before exponentiation, including extreme but finite
        // parameters. This avoids overflowing expm1 or a later int conversion.
        const long double remaining = static_cast<long double>(cap) - base;
        const long double ratio = remaining / scale;
        const long double limit = isfinite(ratio) ? log1pl(ratio)
            : logl(remaining) - logl(static_cast<long double>(scale));
        const long double exponent = static_cast<long double>(rate) * steps;
        if (exponent >= limit) return cap;
        value += static_cast<long double>(scale) * expm1l(exponent);
    }
    if (value >= cap) return cap;
    // Match the fitting convention at floating-point integer boundaries.
    return max(1, static_cast<int>(ceill(value - 1e-10L)));
}
