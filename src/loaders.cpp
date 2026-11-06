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

namespace {
// Two offsets per block keep a human chromosome's index near 1 MB.
const int INDEX_BLOCK = 4096;
const int CURSOR_BLOCK = 1 << 20;

int emission_position(const string &line, const string &path) {
    char *end = nullptr;
    const long position = strtol(line.c_str(), &end, 10);
    if (end == line.c_str() || position < 0 || position > numeric_limits<int>::max())
        throw runtime_error("Invalid or duplicate emission position in " + path);
    return static_cast<int>(position);
}

void parse_emissions(const string &line, const string &path, array<float, NUM_STATES> &row) {
    const char *cursor = line.c_str();
    char *end = nullptr;
    strtol(cursor, &end, 10);
    cursor = end;
    for (int state = 0; state < 5; ++state) {
        row[state] = strtof(cursor, &end);
        if (end == cursor || !isfinite(row[state]))
            throw runtime_error("Expected five finite emission scores in " + path);
        cursor = end;
    }
    row[5] = row[4];
    row[6] = row[4];
    // UniAnn's Perl producer appends the nucleotide for readability.
    // Numerical tables without that optional column are also accepted.
    istringstream rest(cursor);
    string base;
    if (rest >> base) {
        if (base.size() != 1 || string("ACGTRYSWKMBDHVNacgtryswkmbdhvn").find(base[0]) == string::npos)
            throw runtime_error("Invalid optional emission nucleotide in " + path);
        require_end(rest, path);
    }
}

// Record the byte range holding each block's rows and check that every position
// occurs exactly once. Rows may be in any order; a block read scans its whole range,
// which stays short for ascending or descending tables. Scores are validated
// when first loaded.
vector<ScoreSource::ByteRange> index_emissions(const string &path, int length) {
    auto input = open_input(path);
    vector<ScoreSource::ByteRange> ranges(length / INDEX_BLOCK + 1);
    vector<bool> seen(length, false);
    string line;
    uint64_t offset = 0;
    int count = 0;
    while (getline(input, line)) {
        const uint64_t row_offset = offset;
        offset += line.size() + 1;
        if (ignored(line)) continue;
        const int position = emission_position(line, path);
        if (position >= length || seen[position])
            throw runtime_error("Invalid or duplicate emission position in " + path);
        seen[position] = true;
        ++count;
        auto &range = ranges[position / INDEX_BLOCK];
        range.begin = min(range.begin, row_offset);
        range.end = max(range.end, offset);
    }
    if (!count) throw runtime_error("Empty emission file: " + path);
    // A complete emission row is required for every FASTA position.
    if (count < length) {
        const int missing = static_cast<int>(find(seen.begin(), seen.end(), false) - seen.begin());
        throw runtime_error("Missing emission at position " + to_string(missing) + " in " + path);
    }
    return ranges;
}
}

ScoreSource::ScoreSource(const vector<char> &seq, const vector<vector<double>> &trans,
                         const string &emissions, const array<string, 4> &sparse)
    : seq_(seq), trans_(trans), emissions_path_(emissions),
      emission_ranges_(index_emissions(emissions, seq.size())) {
    const int sequence_length = length();
    for (size_t table = 0; table < sparse.size(); ++table) {
        auto input = open_input(sparse[table]);
        vector<pair<int, double>> rows;
        string line;
        while (getline(input, line)) {
            if (ignored(line)) continue;
            istringstream row(line);
            int position;
            double score;
            if (!(row >> position >> score) || position < 0 || position >= sequence_length ||
                !isfinite(score))
                throw runtime_error("Invalid or duplicate motif score in " + sparse[table]);
            require_end(row, sparse[table]);
            rows.emplace_back(position, score);
        }
        stable_sort(rows.begin(), rows.end(),
                    [](const auto &a, const auto &b) { return a.first < b.first; });
        auto &scores = sparse_[table];
        scores.positions.reserve(rows.size());
        scores.values.reserve(rows.size());
        for (const auto &[position, score] : rows) {
            if (!scores.positions.empty() && scores.positions.back() == position)
                throw runtime_error("Invalid or duplicate motif score in " + sparse[table]);
            scores.positions.push_back(position);
            scores.values.push_back(score);
        }
    }
}

ScoreWindow ScoreSource::window(int first, int last) const {
    first = max(first, 0);
    last = min(last, length() - 1);
    ScoreWindow result;
    result.seq_ = &seq_;
    result.trans_ = &trans_;
    result.first_ = first;
    const int count = last - first + 1;
    result.emit_.resize(count);

    auto input = open_input(emissions_path_);
    vector<bool> loaded(count, false);
    int remaining = count;
    string line;
    for (int block = first / INDEX_BLOCK; block <= last / INDEX_BLOCK; ++block) {
        const auto &range = emission_ranges_[block];
        input.clear();
        input.seekg(static_cast<streamoff>(range.begin));
        uint64_t offset = range.begin;
        while (offset < range.end && getline(input, line)) {
            offset += line.size() + 1;
            if (ignored(line)) continue;
            const int position = emission_position(line, emissions_path_);
            if (position < first || position > last || loaded[position - first]) continue;
            parse_emissions(line, emissions_path_, result.emit_[position - first]);
            loaded[position - first] = true;
            --remaining;
        }
    }
    if (remaining) throw runtime_error("Emission file changed while reading: " + emissions_path_);

    for (size_t table = 0; table < sparse_.size(); ++table) {
        const auto &scores = sparse_[table];
        auto &dense = result.sparse_[table];
        dense.assign(count, NEG_INF);
        const auto begin = lower_bound(scores.positions.begin(), scores.positions.end(), first);
        for (size_t i = begin - scores.positions.begin();
             i < scores.positions.size() && scores.positions[i] <= last; ++i)
            dense[scores.positions[i] - first] = scores.values[i];
    }
    return result;
}

ModelInputs ScoreWindow::inputs() const {
    const int count = static_cast<int>(emit_.size());
    return ModelInputs{{emit_.data(), first_, count},
                       {sparse_[0].data(), first_, count},
                       {sparse_[1].data(), first_, count},
                       {sparse_[2].data(), first_, count},
                       {sparse_[3].data(), first_, count},
                       *seq_, *trans_};
}

ModelInputs ScoreCursor::at(int position) {
    if (position + 1 > loaded_last_ && loaded_last_ < source_.length() - 1) {
        window_ = source_.window(position - 2, position + CURSOR_BLOCK);
        loaded_last_ = min(position + CURSOR_BLOCK, source_.length() - 1);
    }
    return window_.inputs();
}
