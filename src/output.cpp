// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#include "altann.hpp"
#include <filesystem>

// Conversion follows motif coordinates, not the raw DP switch position.
// Internal intervals are half-open; GFF features use one-based closed bounds.
// The raw source column remains UniAnn for snapshot compatibility. The CLI's
// export layer maps strand/segment coordinates and adds AltAnn score attributes.
namespace {
// Only state changes affect the annotation, so a run-length path needs no expansion.
class TranscriptBuilder {
public:
    explicit TranscriptBuilder(const string &seqid) : seqid_(seqid) {}

    void step(int from, int to, int dp_position) {
        if (!inside_ && from == 0 && is_exon(to)) {
            inside_ = true;
            current_ = TranscriptAnnotation{};
            current_.sequence_id = seqid_;
            current_.strand = '+';
            current_.path_start = dp_position;
            current_.genomic_start_0based = codon_start_0based(dp_position);
            exon_start_ = current_.genomic_start_0based;
            pending_donor_ = -1;
        }
        if (!inside_) return;

        if (is_exon(from) && is_intron(to)) {
            const int donor = splice_motif_start_0based(dp_position);
            if (exon_start_ < 0 || donor <= exon_start_) {
                inside_ = false;
                return;
            }
            current_.exons.push_back({exon_start_, donor});
            current_.cds_intervals.push_back({exon_start_, donor});
            pending_donor_ = donor;
        } else if (is_intron(from) && is_exon(to)) {
            const int acceptor = splice_motif_start_0based(dp_position);
            if (pending_donor_ < 0 || acceptor < pending_donor_) {
                inside_ = false;
                return;
            }
            current_.junctions.push_back({pending_donor_, acceptor, '+'});
            exon_start_ = acceptor + 2;
            pending_donor_ = -1;
        } else if (is_exon(from) && to == 0) {
            const int transcript_end = dp_position + 1;
            if (exon_start_ < 0 || transcript_end <= exon_start_ ||
                pending_donor_ >= 0) {
                inside_ = false;
                return;
            }
            current_.exons.push_back({exon_start_, transcript_end});
            current_.cds_intervals.push_back({exon_start_, transcript_end});
            current_.path_end = dp_position;
            current_.genomic_end_0based_exclusive = transcript_end;
            transcripts_.push_back(current_);
            inside_ = false;
            exon_start_ = -1;
        }
    }

    vector<TranscriptAnnotation> transcripts_;

private:
    string seqid_;
    bool inside_ = false;
    int exon_start_ = -1;
    int pending_donor_ = -1;
    TranscriptAnnotation current_;
};
}

vector<TranscriptAnnotation> build_transcript_annotations(
    const vector<int> &positions,
    const vector<int> &states,
    const string &seqid
) {
    TranscriptBuilder builder(seqid);
    if (positions.size() != states.size() || states.size() < 2) return builder.transcripts_;
    for (size_t k = 1; k < states.size(); ++k)
        builder.step(states[k - 1], states[k], positions[k]);
    return builder.transcripts_;
}

vector<TranscriptAnnotation> build_transcript_annotations(
    const StatePath &path,
    const string &seqid
) {
    TranscriptBuilder builder(seqid);
    const auto &starts = path.run_starts();
    const auto &states = path.run_states();
    for (size_t run = 1; run < starts.size(); ++run)
        builder.step(states[run - 1], states[run], starts[run]);
    return builder.transcripts_;
}

// Same features and transcript numbering as the UniAnn binary's standard output,
// so UniAnn's wrapper can use this decoder unchanged. Feature scores are written
// as 0; the wrapper's gffread steps discard them.
void write_uniann_gff(ostream &out, const StatePath &path, const string &seqid,
                      const string &fasta_path) {
    const string name = filesystem::path(fasta_path).filename().string();
    int transcript_index = 0;
    const auto feature = [&](int state, int start, int end) {
        if (state == 0) {
            ++transcript_index;
            return;
        }
        out << seqid << "\tUniAnn\t" << (is_exon(state) ? "CDS" : "intron") << '\t' << start << '\t'
            << end << "\t0\t+\t.\tParent=" << name << '.' << transcript_index
            << ";state=" << state_name[state] << '\n';
    };
    out << "##gff-version 3\n";
    const auto &starts = path.run_starts();
    const auto &states = path.run_states();
    int previous = states[0], start = 0;
    for (size_t run = 1; run < starts.size(); ++run) {
        const int current = states[run - 1], next = states[run], position = starts[run];
        int end = position - 1;
        if (current == 0) end = position - 2;
        if (next == 0) end = position + 1;
        if (end > start) {
            if (is_exon(current) && is_intron(previous)) feature(current, start + 2, end);
            else if (is_intron(current) && is_exon(previous)) feature(current, start, end + 2);
            else feature(current, start, end);
        }
        if (current == 0) start = position - 1;
        else if (next == 0) start = position + 2;
        else start = position;
        previous = current;
    }
    feature(states.back(), start, path.size() - 1);
}

string intron_chain_key(const TranscriptAnnotation &transcript) {
    ostringstream out;
    out << transcript.sequence_id << '|' << transcript.strand << '|';
    for (const auto &junction : transcript.junctions)
        out << junction.donor_0based << '-' << junction.acceptor_0based << ',';
    return out.str();
}

string full_transcript_key(const TranscriptAnnotation &transcript) {
    ostringstream out;
    out << intron_chain_key(transcript) << '|' << transcript.genomic_start_0based
        << '-' << transcript.genomic_end_0based_exclusive << '|';
    for (const auto &exon : transcript.exons)
        out << exon.start_0based << '-' << exon.end_0based_exclusive << ',';
    out << '|';
    for (const auto &cds : transcript.cds_intervals)
        out << cds.start_0based << '-' << cds.end_0based_exclusive << ',';
    return out.str();
}

static bool interval_contains(const GenomicInterval &outer,
                              const GenomicInterval &inner) {
    return outer.start_0based <= inner.start_0based &&
           inner.end_0based_exclusive <= outer.end_0based_exclusive;
}

string classify_splice_difference(
    const TranscriptAnnotation &reference,
    const TranscriptAnnotation &candidate
) {
    if (reference.junctions == candidate.junctions) return "identical";
    if (reference.junctions.size() == candidate.junctions.size()) {
        int donor_changes = 0;
        int acceptor_changes = 0;
        for (size_t i = 0; i < reference.junctions.size(); ++i) {
            donor_changes += reference.junctions[i].donor_0based !=
                             candidate.junctions[i].donor_0based;
            acceptor_changes += reference.junctions[i].acceptor_0based !=
                                candidate.junctions[i].acceptor_0based;
        }
        if (donor_changes == 1 && acceptor_changes == 0)
            return "alternative_donor";
        if (donor_changes == 0 && acceptor_changes == 1)
            return "alternative_acceptor";
    }

    for (size_t i = 1; i + 1 < reference.exons.size(); ++i) {
        bool covered = false;
        for (const auto &exon : candidate.exons)
            covered = covered || interval_contains(exon, reference.exons[i]);
        if (!covered && candidate.exons.size() < reference.exons.size())
            return "exon_skipping";
    }
    for (size_t i = 1; i + 1 < candidate.exons.size(); ++i) {
        bool covered = false;
        for (const auto &exon : reference.exons)
            covered = covered || interval_contains(exon, candidate.exons[i]);
        if (!covered && candidate.exons.size() > reference.exons.size())
            return "alternative_exon";
    }
    for (const auto &junction : reference.junctions) {
        GenomicInterval intron{junction.donor_0based,
                               junction.acceptor_0based + 2};
        for (const auto &exon : candidate.exons)
            if (interval_contains(exon, intron)) return "intron_retention";
    }
    return "complex";
}

// Donor and acceptor positions are compared as separate labelled sets, so a
// shifted splice site breaks only its own boundary rather than whole introns.
// Two intronless models are identical under this metric.
double splice_site_jaccard(const TranscriptAnnotation &reference,
                           const TranscriptAnnotation &candidate) {
    if (reference.sequence_id != candidate.sequence_id ||
        reference.strand != candidate.strand) return 0.0;
    auto sites = [](const TranscriptAnnotation &transcript, bool donor) {
        vector<int> result;
        for (const auto &junction : transcript.junctions)
            result.push_back(donor ? junction.donor_0based : junction.acceptor_0based);
        sort(result.begin(), result.end());
        result.erase(unique(result.begin(), result.end()), result.end());
        return result;
    };
    size_t shared = 0, total = 0;
    for (bool donor : {true, false}) {
        const auto a = sites(reference, donor), b = sites(candidate, donor);
        vector<int> common;
        set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(common));
        shared += common.size();
        total += a.size() + b.size() - common.size();
    }
    return total == 0 ? 1.0 : static_cast<double>(shared) / total;
}

static string format_score(double score) {
    if (!isfinite(score) || score <= NEG_INF) return ".";
    ostringstream out;
    out << setprecision(17) << score;
    return out.str();
}

void write_hierarchical_transcript(
    ostream &out,
    const TranscriptAnnotation &transcript,
    const string &gene_id,
    const string &transcript_id
) {
    out << transcript.sequence_id << "\tUniAnn\ttranscript\t"
        << transcript.genomic_start_0based + 1 << '\t'
        << transcript.genomic_end_0based_exclusive << '\t'
        << format_score(transcript.score) << "\t+\t.\tID=" << transcript_id
        << ";Parent=" << gene_id << "\n";

    int cumulative_cds = 0;
    for (size_t i = 0; i < transcript.exons.size(); ++i) {
        const auto &exon = transcript.exons[i];
        out << transcript.sequence_id << "\tUniAnn\texon\t"
            << exon.start_0based + 1 << '\t' << exon.end_0based_exclusive
            << "\t.\t+\t.\tID=" << transcript_id << ".exon" << i + 1
            << ";Parent=" << transcript_id << "\n";
        const auto &cds = transcript.cds_intervals[i];
        const int phase = i == 0 ? 0 : (3 - cumulative_cds % 3) % 3;
        out << transcript.sequence_id << "\tUniAnn\tCDS\t"
            << cds.start_0based + 1 << '\t' << cds.end_0based_exclusive
            << "\t.\t+\t" << phase << "\tID=" << transcript_id
            << ".cds" << i + 1 << ";Parent=" << transcript_id << "\n";
        cumulative_cds += cds.end_0based_exclusive - cds.start_0based;
        if (i < transcript.junctions.size()) {
            const auto &junction = transcript.junctions[i];
            out << transcript.sequence_id << "\tUniAnn\tintron\t"
                << junction.donor_0based + 1 << '\t'
                << junction.acceptor_0based + 2
                << "\t.\t+\t.\tID=" << transcript_id << ".intron"
                << i + 1 << ";Parent=" << transcript_id << "\n";
        }
    }
}

string tsv_escape_chain(const TranscriptAnnotation &transcript) {
    ostringstream out;
    for (size_t i = 0; i < transcript.junctions.size(); ++i) {
        if (i) out << ',';
        out << transcript.junctions[i].donor_0based << '-'
            << transcript.junctions[i].acceptor_0based;
    }
    return out.str();
}
