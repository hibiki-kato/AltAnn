// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#include "altann.hpp"

// Transition rules preserve the Dmel seven-state decoder, including its finite
// penalties and motif offsets. All model coordinates are zero-based.
vector<vector<double>> init_transitions() {
    vector<vector<double>> trans(NUM_STATES, vector<double>(NUM_STATES, NEG_INF));

    double self_prob = 0.0;

    // Noncoding self
    trans[0][0] = self_prob;

    // Exon frame cycling: none (self only)
    trans[1][1] = self_prob;
    trans[2][2] = self_prob;
    trans[3][3] = self_prob;

    // Intron states: self only
    trans[4][4] = self_prob;
    trans[5][5] = self_prob;
    trans[6][6] = self_prob;

    return trans;
}

int splice_motif_start_0based(int dp_position) {
    return dp_position - 1;
}

int codon_start_0based(int dp_position) {
    return std::max(0, dp_position - 2);
}

static bool sequence_has_dinucleotide(const vector<char> &seq, int start,
                                      char first, char second) {
    return start >= 0 && start + 1 < static_cast<int>(seq.size()) &&
           toupper(seq[start]) == first && toupper(seq[start + 1]) == second;
}

static string sequence_codon_ending_at(const vector<char> &seq, int position) {
    if (position < 2 || position >= static_cast<int>(seq.size())) return "";
    string codon;
    codon.push_back(toupper(seq[position - 2]));
    codon.push_back(toupper(seq[position - 1]));
    codon.push_back(toupper(seq[position]));
    return codon;
}

static bool is_stop_codon(const string &codon) {
    return codon == "TAA" || codon == "TAG" || codon == "TGA";
}

static double destination_emission(int position, int to,
                                   const ModelInputs &inputs) {
    double emission = inputs.emit[position][to];
    // Preserve the legacy TAG/AG overlap convention exactly.
    if (emission <= -1e6 && position > 0 && position + 1 <
            static_cast<int>(inputs.seq.size()) &&
        toupper(inputs.seq[position - 1]) == 'A' &&
        toupper(inputs.seq[position]) == 'G') {
        emission = inputs.emit[position + 1][to];
    }
    return emission;
}

TransitionResult evaluate_transition(
    int position,
    int from,
    int to,
    const PathMetadata &previous,
    const ModelInputs &inputs
) {
    TransitionResult result;
    if (position <= 0 || position >= static_cast<int>(inputs.seq.size()) ||
        from < 0 || from >= NUM_STATES || to < 0 || to >= NUM_STATES) {
        return result;
    }

    double log_t = inputs.trans[from][to];
    const string codon = sequence_codon_ending_at(inputs.seq, position);
    const bool stop = is_stop_codon(codon);

    if (is_exon(to) && is_exon(from) && from == to && stop &&
        codon_start_0based(position) % 3 == to - 1) {
        log_t = NEG_INF;
    }

    if (is_exon(from) && is_intron(to) &&
        sequence_has_dinucleotide(inputs.seq,
                                  splice_motif_start_0based(position), 'G', 'T')) {
        const int length = previous.exon_len - 2;
        log_t = NEG_INF;
        if ((to - 4) == (from - 1) &&
            (length >= MIN_EXON || previous.exon_from == 0)) {
            log_t = inputs.gt_score[splice_motif_start_0based(position)];
        }
    }

    if (is_intron(from) && is_exon(to) &&
        sequence_has_dinucleotide(inputs.seq,
                                  splice_motif_start_0based(position), 'A', 'G')) {
        const int length = previous.intron_len + 2;
        log_t = NEG_INF;
        if (length >= MIN_INTRON) {
            const int intron_frame = from - 4;
            const int exon_frame = to - 1;
            if ((exon_frame - intron_frame + 3) % 3 == length % 3) {
                log_t = inputs.ag_score[splice_motif_start_0based(position)];
            }
        }
    }

    if (is_exon(from) && to == 0 && stop) {
        const int length = previous.exon_len - 2;
        const int frame = from - 1;
        if (codon_start_0based(position) % 3 == frame) {
            if (is_intron(previous.exon_from) ||
                (previous.exon_from == 0 && length > MIN_SINGLE)) {
                log_t = inputs.stop_score[codon_start_0based(position)];
            } else {
                // This finite legacy penalty is a model rule, not an invalid
                // transition sentinel, and is retained for compatibility.
                log_t = -1e3;
            }
        }
    }

    if (is_exon(to) && from == 0 && position >= 2) {
        const int length = previous.inter_len + 2;
        if (codon == "ATG" &&
            (length >= MIN_INTER || position < MIN_INTER) &&
            previous.predecessor_state == 0) {
            const int frame = to - 1;
            if (codon_start_0based(position) % 3 == frame) {
                log_t = position < 25
                    ? 1.0
                    : inputs.atg_score[codon_start_0based(position)];
            }
        }
    }

    result.transition_score = log_t;
    result.emission_score = destination_emission(position, to, inputs);
    result.allowed = log_t > NEG_INF && result.emission_score > NEG_INF;

    PathMetadata next;
    next.predecessor_state = from;
    if (is_intron(to))
        next.intron_len = is_intron(from) ? previous.intron_len + 1 : 1;
    if (is_exon(to)) {
        if (is_exon(from)) {
            next.exon_len = previous.exon_len + 1;
            next.exon_from = previous.exon_from;
        } else {
            next.exon_len = from == 0 ? 3 : 1;
            next.exon_from = from;
        }
    }
    if (to == 0)
        next.inter_len = from == 0 ? previous.inter_len + 1 : 1;
    result.next_metadata = next;
    return result;
}
