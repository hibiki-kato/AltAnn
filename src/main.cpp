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
        << "  [--raw] Write every gene of every decoded path, with no overlap filter or deduplication.\n"
        << "  [--viterbi-log FILE] Read the saved UniAnn dp/bt baseline instead of rerunning it.\n"
        << "  [--score-models FILE] Score supplied CDS models instead of K-best decoding;\n"
        << "    FILE rows are ID<TAB>start-end,... (zero-based half-open); --report is unused.\n"
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
            if (name == "--raw") { options.raw = true; continue; }
            if (arg + 1 >= argc) throw runtime_error("Missing value for " + name);
            const string value = argv[++arg];
            if (name == "--output") options.output_filename = value;
            else if (name == "--report") options.report_filename = value;
            else if (name == "--viterbi-log") options.viterbi_log_filename = value;
            else if (name == "--score-models") options.models_filename = value;
            else if (name == "--k") options.k = integer_option(value, name);
            else if (name == "--flank") options.flank = integer_option(value, name);
            else if (name == "--threads") options.threads = value == "auto" ? 0 : integer_option(value, name);
            else throw runtime_error("Unknown option: " + name);
        }
        if (options.k < 1 || options.k > 32767) throw runtime_error("K must be between 1 and 32767");
        if (options.flank < -1 || options.flank > numeric_limits<int>::max() / 2)
            throw runtime_error("Flank must be -1 or a nonnegative supported distance");
        if (options.threads < 0) throw runtime_error("Threads must be nonnegative");
        const bool scoring = !options.models_filename.empty();
        if (options.output_filename.empty() || (!scoring && options.report_filename.empty()))
            throw runtime_error(scoring ? "--output is required"
                                        : "Both --output and --report are required");
        if (scoring && options.raw)
            throw runtime_error("--raw cannot be combined with --score-models");
        if (scoring && !options.report_filename.empty())
            throw runtime_error("--report cannot be combined with --score-models");
        // Resolve existing symlinks and relative components before opening
        // outputs, so a typo cannot replace an input or alias another output.
        set<filesystem::path> output_paths;
        vector<string> input_paths(argv + 1, argv + 7);
        if (!options.viterbi_log_filename.empty())
            input_paths.push_back(options.viterbi_log_filename);
        if (scoring) input_paths.push_back(options.models_filename);
        vector<string> outputs{options.output_filename};
        if (!scoring) {
            outputs.push_back(options.report_filename);
            outputs.push_back(options.output_filename + ".reference.gff");
        }
        for (const string &path : outputs) {
            const auto resolved = filesystem::weakly_canonical(filesystem::absolute(path));
            if (!output_paths.insert(resolved).second)
                throw runtime_error("Output files must have distinct paths");
            for (const string &input_path : input_paths) {
                const auto input = filesystem::weakly_canonical(filesystem::absolute(input_path));
                if (resolved == input || (filesystem::exists(path) && filesystem::exists(input)
                                         && filesystem::equivalent(path, input)))
                    throw runtime_error("Output path aliases an input: " + path);
            }
        }
        const auto fasta = read_fasta(argv[1]);
        const auto transitions = init_transitions();
        const ScoreSource source(fasta.sequence, transitions, argv[2],
                                 {argv[3], argv[4], argv[5], argv[6]});
        const auto baseline = options.viterbi_log_filename.empty()
            ? decode_baseline(source)
            : load_baseline_log(source, options.viterbi_log_filename);
        const auto references = build_transcript_annotations(baseline.states, fasta.id);
        if (scoring) {
            run_model_scoring(fasta.id, references, baseline, source, options);
            cerr << "Scored models against " << references.size() << " reference loci\n";
            return 0;
        }
        if (!run_local_k_best(fasta.id, references, baseline, source, options)) return 1;
        cerr << "Decoded " << references.size() << " reference loci\n";
        return 0;
    } catch (const exception &error) {
        cerr << "altann-core: " << error.what() << '\n';
        return 1;
    }
}
