#include "altann.hpp"
#include <filesystem>

namespace {
int integer_option(const string &value, const string &name) {
    size_t end = 0;
    int number;
    try { number = stoi(value, &end); }
    catch (const exception &) { throw runtime_error("Invalid integer for " + name); }
    if (end != value.size()) throw runtime_error("Invalid integer for " + name);
    return number;
}
void usage(ostream &out) {
    out << "Usage: altann-core FASTA EMISSIONS GT AG ATG STOP\n"
        << "  --output FILE --report FILE [--k 10] [--flank 1000] [--threads 0]\n"
        << "Coordinates and sequence orientation follow the input FASTA.\n"
        << "Threads 0 selects available CPUs; flank -1 selects whole intergenic bounds.\n";
}
}
int main(int argc, char **argv) {
    try {
        if (argc == 2 && string(argv[1]) == "--help") { usage(cout); return 0; }
        if (argc < 7) { usage(cerr); return 2; }
        Options options;
        for (int arg = 7; arg < argc; ++arg) {
            const string name = argv[arg];
            if (arg + 1 >= argc) throw runtime_error("Missing value for " + name);
            const string value = argv[++arg];
            if (name == "--output") options.output_filename = value;
            else if (name == "--report") options.report_filename = value;
            else if (name == "--k") options.k = integer_option(value, name);
            else if (name == "--flank") options.flank = integer_option(value, name);
            else if (name == "--threads") options.threads = value == "auto" ? 0 : integer_option(value, name);
            else throw runtime_error("Unknown option: " + name);
        }
        if (options.k < 1 || options.k > 32767) throw runtime_error("K must be between 1 and 32767");
        if (options.flank < -1 || options.flank > numeric_limits<int>::max() / 2)
            throw runtime_error("Flank must be -1 or a nonnegative supported distance");
        if (options.threads < 0) throw runtime_error("Threads must be nonnegative");
        if (options.output_filename.empty() || options.report_filename.empty())
            throw runtime_error("Both --output and --report are required");
        // Resolve existing symlinks and relative components before opening
        // outputs, so a typo cannot replace an input or alias another output.
        set<filesystem::path> output_paths;
        for (const string &path : {options.output_filename, options.report_filename,
                                  options.output_filename + ".reference.gff"}) {
            const auto resolved = filesystem::weakly_canonical(filesystem::absolute(path));
            if (!output_paths.insert(resolved).second)
                throw runtime_error("Output files must have distinct paths");
            for (int arg = 1; arg <= 6; ++arg) {
                const auto input = filesystem::weakly_canonical(filesystem::absolute(argv[arg]));
                if (resolved == input || (filesystem::exists(path) && filesystem::exists(input)
                                         && filesystem::equivalent(path, input)))
                    throw runtime_error("Output path aliases an input: " + path);
            }
        }
        const auto fasta = read_fasta(argv[1]);
        const int length = fasta.sequence.size();
        const auto emissions = load_emissions(argv[2], length);
        const auto gt = load_sparse_scores(argv[3], length);
        const auto ag = load_sparse_scores(argv[4], length);
        const auto atg = load_sparse_scores(argv[5], length);
        const auto stop = load_sparse_scores(argv[6], length);
        const auto transitions = init_transitions();
        const ModelInputs inputs{emissions, gt, ag, atg, stop, fasta.sequence, transitions};
        const auto baseline = decode_baseline(inputs);
        const auto references = build_transcript_annotations(baseline.states, fasta.id);
        if (!run_local_k_best(fasta.id, references, baseline.states,
                              baseline.noncoding_prefix, inputs, options)) return 1;
        cerr << "Decoded " << references.size() << " reference loci\n";
        return 0;
    } catch (const exception &error) {
        cerr << "altann-core: " << error.what() << '\n';
        return 1;
    }
}
