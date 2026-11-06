// Exercise model-to-state conversion independently of input files.
#include "../src/score.cpp"

namespace {
void require(bool condition, const string &message) {
    if (!condition) throw runtime_error(message);
}

vector<GenomicInterval> round_trip(const vector<GenomicInterval> &cds, int length) {
    int first = -1;
    vector<int> states;
    require(model_states(cds, length, first, states).empty(), "model converts");
    vector<int> path(length, 0);
    for (size_t i = 0; i < states.size(); ++i) path[first + i] = states[i];
    vector<int> positions(path.size());
    for (size_t i = 0; i < positions.size(); ++i) positions[i] = i;
    const auto annotations = build_transcript_annotations(positions, path, "fixture");
    require(annotations.size() == 1, "one transcript is rebuilt");
    return annotations[0].cds_intervals;
}

void test_round_trip_and_frames() {
    const vector<GenomicInterval> single{{100, 400}};
    require(round_trip(single, 1000) == single, "single exon round trip");
    const vector<GenomicInterval> spliced{{101, 200}, {300, 451}, {600, 749}};
    require(round_trip(spliced, 1000) == spliced, "spliced round trip");
    int first = -1;
    vector<int> states;
    model_states(spliced, 1000, first, states);
    require(first == 103, "N->E switch two bases after the CDS start");
    require(states[0] == 1 + 101 % 3, "first exon frame is the start codon phase");
    require(states[200 - first] == states[0], "donor G is the last exon state");
    require(states[201 - first] == 4 + 101 % 3, "intron keeps the preceding exon frame");
    const int frame2 = (101 % 3 + 300 - 200) % 3;
    require(states[299 - first] == 1 + frame2, "acceptor G starts the next exon frame");
    require(states.back() == states[600 - 1 - first], "last exon keeps one frame");
    require(first + static_cast<int>(states.size()) - 1 == 747, "stop codon's last base is N");
}

void test_rejected_structures() {
    int first = -1;
    vector<int> states;
    require(model_states({{100, 200}, {201, 300}}, 1000, first, states) ==
            "overlapping_or_adjacent_cds", "adjacent CDS is rejected");
    require(model_states({{100, 103}}, 1000, first, states) == "cds_too_short", "tiny CDS is rejected");
    require(model_states({{700, 1000}}, 1000, first, states).empty(),
            "a stop ending at the last base still has an N anchor");
    require(model_states({{700, 1001}}, 1000, first, states) == "touches_sequence_end",
            "CDS beyond the sequence is rejected");
}
} // namespace

int main() {
    try {
        test_round_trip_and_frames();
        test_rejected_structures();
    } catch (const exception &error) {
        cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
