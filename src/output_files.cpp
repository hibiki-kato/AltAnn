#include "altann.hpp"
#include <filesystem>
#include <cstdlib>

namespace fs = std::filesystem;
namespace {
struct StagedFile {
    fs::path target, directory;
    bool previous_saved = false;
    bool installed = false;
};
}

void write_output_bundle(const vector<pair<string, string>> &files) {
    vector<StagedFile> staged;
    staged.reserve(files.size());
    try {
        // Stage on the target filesystem so each final rename is atomic. No
        // destination changes until every file has been fully written/closed.
        for (const auto &file : files) {
            StagedFile entry;
            entry.target = fs::absolute(file.first);
            if (fs::is_directory(entry.target))
                throw runtime_error("Output path is a directory: " + file.first);
            string pattern = entry.target.string() + ".altann-XXXXXX";
            vector<char> buffer(pattern.begin(), pattern.end());
            buffer.push_back('\0');
            const char *directory = mkdtemp(buffer.data());
            if (!directory) throw runtime_error("Cannot stage output: " + file.first);
            entry.directory = directory;
            staged.push_back(move(entry));
            ofstream output(staged.back().directory / "new", ios::binary);
            output.exceptions(ios::badbit | ios::failbit);
            output.write(file.second.data(), file.second.size());
            output.close();
        }
        for (auto &entry : staged) {
            if (fs::exists(entry.target) || fs::is_symlink(entry.target)) {
                fs::rename(entry.target, entry.directory / "previous");
                entry.previous_saved = true;
            }
            fs::rename(entry.directory / "new", entry.target);
            entry.installed = true;
        }
    } catch (...) {
        // Restore prior outputs after ordinary write/rename failures. Retain a
        // backup directory if restoration itself fails (e.g. filesystem loss).
        for (auto it = staged.rbegin(); it != staged.rend(); ++it) {
            error_code error;
            bool restored = true;
            if (it->installed) {
                fs::remove(it->target, error);
                restored = !error;
            }
            if (it->previous_saved && restored) {
                fs::rename(it->directory / "previous", it->target, error);
                restored = !error;
            }
            if (restored) fs::remove_all(it->directory, error);
            else cerr << "Output recovery backup retained at " << it->directory << '\n';
        }
        throw;
    }
    for (const auto &entry : staged) {
        error_code error;
        fs::remove_all(entry.directory, error);
        if (error) cerr << "Cannot remove staging directory " << entry.directory << '\n';
    }
}
