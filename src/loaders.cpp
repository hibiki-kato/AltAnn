#include "altann.hpp"

// Read one oriented sequence and its zero-based numerical score tables.
// Motif tables are sparse; emissions must cover every sequence position.

namespace {
ifstream open_input(const string &path) {
    ifstream input(path);
    if (!input) throw runtime_error("Cannot open input: " + path);
    return input;
}
bool ignored(const string &line) {
    const auto start = line.find_first_not_of(" \t\r\n");
    return start == string::npos || line[start] == '#';
}
void require_end(istringstream &row, const string &path) {
    string extra;
    if (row >> extra) throw runtime_error("Unexpected extra column in " + path);
}
}

Fasta read_fasta(const string &path) {
    auto input = open_input(path);
    Fasta fasta;
    string line;
    bool header_seen = false;
    while (getline(input, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (header_seen) throw runtime_error("Expected a single FASTA record: " + path);
            header_seen = true;
            istringstream header(line.substr(1));
            if (!(header >> fasta.id)) throw runtime_error("Empty FASTA identifier: " + path);
            continue;
        }
        for (unsigned char base : line) {
            if (isspace(base)) continue;
            if (!header_seen) throw runtime_error("Sequence before FASTA header: " + path);
            base = toupper(base);
            if (string("ACGTRYSWKMBDHVN").find(base) == string::npos)
                throw runtime_error("Invalid DNA symbol in " + path);
            if (fasta.sequence.size() >= static_cast<size_t>(numeric_limits<int>::max() - 2))
                throw runtime_error("FASTA record exceeds supported coordinate range");
            fasta.sequence.push_back(base);
        }
    }
    if (fasta.sequence.empty()) throw runtime_error("Empty FASTA sequence: " + path);
    return fasta;
}

vector<array<float, NUM_STATES>> load_emissions(const string &path, int length) {
    auto input = open_input(path);
    vector<array<float, NUM_STATES>> emissions(length);
    for (auto &row : emissions) row.fill(NEG_INF);
    vector<bool> seen(length, false);
    string line;
    int count = 0;
    while (getline(input, line)) {
        if (ignored(line)) continue;
        istringstream row(line);
        int position;
        if (!(row >> position) || position < 0 || position >= length || seen[position])
            throw runtime_error("Invalid or duplicate emission position in " + path);
        for (int state = 0; state < 5; ++state) {
            if (!(row >> emissions[position][state]) || !isfinite(emissions[position][state]))
                throw runtime_error("Expected five finite emission scores in " + path);
        }
        // UniAnn's Perl producer appends the nucleotide for readability.
        // Numerical tables without that optional column are also accepted.
        string base;
        if (row >> base) {
            if (base.size() != 1 || string("ACGTRYSWKMBDHVNacgtryswkmbdhvn").find(base[0]) == string::npos)
                throw runtime_error("Invalid optional emission nucleotide in " + path);
            require_end(row, path);
        }
        emissions[position][5] = emissions[position][4];
        emissions[position][6] = emissions[position][4];
        seen[position] = true;
        ++count;
    }
    if (!count) throw runtime_error("Empty emission file: " + path);
    // A complete emission row is required for every FASTA position.
    for (int position = 0; position < length; ++position)
        if (!seen[position]) throw runtime_error("Missing emission at position " + to_string(position) + " in " + path);
    return emissions;
}

vector<double> load_sparse_scores(const string &path, int length) {
    auto input = open_input(path);
    vector<double> scores(length, NEG_INF);
    vector<bool> seen(length, false);
    string line;
    while (getline(input, line)) {
        if (ignored(line)) continue;
        istringstream row(line);
        int position;
        double score;
        if (!(row >> position >> score) || position < 0 || position >= length ||
            !isfinite(score) || seen[position])
            throw runtime_error("Invalid or duplicate motif score in " + path);
        require_end(row, path);
        scores[position] = score;
        seen[position] = true;
    }
    return scores;
}
