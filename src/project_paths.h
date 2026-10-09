#ifndef PROJECT_PATHS_H
#define PROJECT_PATHS_H

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#include <cstdlib>
#include <vector>
#elif defined(__linux__)
#include <limits.h>
#include <unistd.h>
#endif

// All paths returned by this header are absolute. The project root is found
// from the executable, source file, or current directory, so the project does
// not depend on a particular Windows username or installation directory.
namespace project_paths {
namespace fs = std::filesystem;

inline fs::path executable_path() {
#ifdef _WIN32
    std::vector<wchar_t> buffer(32768, L'\0');
    const DWORD length = GetModuleFileNameW(nullptr, buffer.data(),
                                             static_cast<DWORD>(buffer.size()));
    if (length > 0 && length < buffer.size()) {
        return fs::path(std::wstring(buffer.data(), length));
    }
#elif defined(__APPLE__)
    uint32_t size = 0;
    if (_NSGetExecutablePath(nullptr, &size) == -1 && size > 0) {
        std::vector<char> buffer(size + 1, '\0');
        if (_NSGetExecutablePath(buffer.data(), &size) == 0) {
            std::error_code ec;
            fs::path p(buffer.data());
            fs::path canonical = fs::weakly_canonical(p, ec);
            return ec ? fs::absolute(p) : canonical;
        }
    }
#elif defined(__linux__)
    std::vector<char> buffer(PATH_MAX + 1, '\0');
    const ssize_t length = readlink("/proc/self/exe", buffer.data(), PATH_MAX);
    if (length > 0) {
        buffer[static_cast<std::size_t>(length)] = '\0';
        return fs::path(buffer.data());
    }
#endif
    return {};
}

inline bool is_project_root(const fs::path& p) {
    std::error_code ec;
    return fs::is_directory(p / "data", ec) &&
           fs::is_directory(p / "src", ec) &&
           fs::is_directory(p / "scripts", ec);
}

inline fs::path search_upwards(fs::path start) {
    if (start.empty()) return {};

    std::error_code ec;
    start = fs::absolute(start, ec);
    if (ec) return {};

    if (!fs::is_directory(start, ec)) {
        start = start.parent_path();
    }

    while (!start.empty()) {
        if (is_project_root(start)) {
            fs::path resolved = fs::weakly_canonical(start, ec);
            return ec ? start : resolved;
        }

        const fs::path parent = start.parent_path();
        if (parent == start) break;
        start = parent;
    }
    return {};
}

inline const fs::path& root() {
    static const fs::path project_root = []() -> fs::path {
        std::vector<fs::path> candidates;
        const fs::path exe = executable_path();
        if (!exe.empty()) candidates.push_back(exe);
        candidates.emplace_back(fs::path(__FILE__));
        candidates.push_back(fs::current_path());

        for (const auto& candidate : candidates) {
            const fs::path found = search_upwards(candidate);
            if (!found.empty()) return found;
        }

        throw std::runtime_error(
            "Could not locate the project root. Expected folders: data, src, and scripts.");
    }();
    return project_root;
}

inline fs::path data_file(const std::string& filename) {
    return root() / "data" / filename;
}

inline fs::path script_file(const std::string& filename) {
    return root() / "scripts" / filename;
}

inline fs::path output_dir() {
    const fs::path dir = root() / "outputs";
    std::error_code ec;
    fs::create_directories(dir, ec);
    if (ec) {
        throw std::runtime_error("Could not create output directory: " + dir.string() +
                                 " (" + ec.message() + ")");
    }
    return dir;
}

inline fs::path output_file(const std::string& filename) {
    return output_dir() / filename;
}

} // namespace project_paths

#endif // PROJECT_PATHS_H
