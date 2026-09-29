// Adapted from UniAnn by Aleksey Zimin and contributors. See THIRD_PARTY.md.
#include "altann.hpp"

// Conversion follows motif coordinates, not the raw DP switch position.
// Internal intervals are half-open; GFF features use one-based closed bounds.
// The raw source column remains UniAnn for snapshot compatibility. The CLI's
// export layer maps strand/segment coordinates and adds AltAnn score attributes.
vector<TranscriptAnnotation> build_transcript_annotations(
    const vector<int> &positions,
    const vector<int> &states,
    const string &seqid
) {
    vector<TranscriptAnnotation> transcripts;
    if (positions.size() != states.size() || states.size() < 2) return transcripts;

    bool inside = false;
    int exon_start = -1;
    int pending_donor = -1;
    TranscriptAnnotation current;
    current.sequence_id = seqid;

    for (size_t k = 1; k < states.size(); ++k) {
        const int from = states[k - 1];
        const int to = states[k];
        const int dp_position = positions[k];

        if (!inside && from == 0 && is_exon(to)) {
            inside = true;
            current = TranscriptAnnotation{};
            current.sequence_id = seqid;
            current.strand = '+';
            current.path_start = dp_position;
            current.genomic_start_0based = codon_start_0based(dp_position);
            exon_start = current.genomic_start_0based;
            pending_donor = -1;
        }
        if (!inside) continue;

        if (is_exon(from) && is_intron(to)) {
            const int donor = splice_motif_start_0based(dp_position);
            if (exon_start < 0 || donor <= exon_start) {
                inside = false;
                continue;
            }
            current.exons.push_back({exon_start, donor});
            current.cds_intervals.push_back({exon_start, donor});
            pending_donor = donor;
        } else if (is_intron(from) && is_exon(to)) {
            const int acceptor = splice_motif_start_0based(dp_position);
            if (pending_donor < 0 || acceptor < pending_donor) {
                inside = false;
                continue;
            }
            current.junctions.push_back({pending_donor, acceptor, '+'});
            exon_start = acceptor + 2;
            pending_donor = -1;
        } else if (is_exon(from) && to == 0) {
            const int transcript_end = dp_position + 1;
            if (exon_start < 0 || transcript_end <= exon_start ||
                pending_donor >= 0) {
                inside = false;
                continue;
            }
            current.exons.push_back({exon_start, transcript_end});
            current.cds_intervals.push_back({exon_start, transcript_end});
            current.path_end = dp_position;
            current.genomic_end_0based_exclusive = transcript_end;
            transcripts.push_back(current);
            inside = false;
            exon_start = -1;
        }
    }
    return transcripts;
}

vector<TranscriptAnnotation> build_transcript_annotations(
    const vector<int> &states,
    const string &seqid
) {
    vector<int> positions(states.size());
    for (size_t i = 0; i < positions.size(); ++i) positions[i] = i;
    return build_transcript_annotations(positions, states, seqid);
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
