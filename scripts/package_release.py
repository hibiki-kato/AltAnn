#!/usr/bin/env python3
"""Bundle an installed AltAnn tree and its OpenMP runtime for distribution.

The build machine needs the native toolchain and runtime development package.
Archive users need Python 3.9 or later. No compiler, Perl, PSAURON, or
separately installed OpenMP runtime is needed.

Run after cmake --install: package_release.py --prefix INSTALL --output FILE.tar.gz
The installed tree and repository are never modified. All binary rewriting
happens in a private staging directory. No dependencies or licenses are fetched.
"""

import argparse
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile


ROOT = Path(__file__).resolve().parents[1]
LINUX_RUNTIMES = ("libgomp.so.1", "libstdc++.so.6", "libgcc_s.so.1")


def run(*arguments, env=None):
    """Run a build-machine inspection tool and retain readable failure context."""
    try:
        result = subprocess.run([str(arg) for arg in arguments], check=True,
                                capture_output=True, text=True, env=env)
    except FileNotFoundError as error:
        raise RuntimeError(f"Required packaging tool is missing: {arguments[0]}") from error
    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or "").strip()
        raise RuntimeError(f"{arguments[0]} failed: {detail}") from error
    return result.stdout


def clean_runtime_environment():
    """Prevent a developer's loader overrides from hiding an incomplete bundle."""
    return {key: value for key, value in os.environ.items()
            if key not in {"LD_LIBRARY_PATH", "LD_PRELOAD", "DYLD_LIBRARY_PATH",
                           "DYLD_FALLBACK_LIBRARY_PATH", "DYLD_INSERT_LIBRARIES"}}


def copy_license(candidates, destination, instruction):
    """Use an installed license, with an explicit failure instead of a download."""
    for candidate in candidates:
        candidate = Path(candidate)
        if candidate.is_file() and candidate.stat().st_size:
            shutil.copy2(candidate, destination)
            return candidate
    raise RuntimeError(f"Runtime license is unavailable. {instruction}")


def linux_dependencies(binary):
    text = run("ldd", binary, env=clean_runtime_environment())
    if "not found" in text:
        raise RuntimeError(f"Unresolved native dependencies for {binary.name}:\n{text}")
    dependencies = {}
    for line in text.splitlines():
        match = re.match(r"\s*(\S+)\s+=>\s+(/.*?)\s+\(0x[0-9a-fA-F]+\)", line)
        if match:
            dependencies[match.group(1)] = Path(match.group(2))
    return dependencies


def bundle_linux(tree, licenses):
    """Bundle GNU runtimes and require the binary's relative RPATH to find them."""
    core = tree / "bin" / "altann-core"
    dependencies = linux_dependencies(core)
    if "libgomp.so.1" not in dependencies:
        raise RuntimeError("Expected GNU OpenMP (libgomp.so.1). Build this Linux release with GCC/OpenMP.")
    lib = tree / "lib"
    lib.mkdir(exist_ok=True)
    bundled = []
    for name in LINUX_RUNTIMES:
        if name in dependencies:
            source = dependencies[name].resolve()
            if source != (lib / name).resolve():
                shutil.copy2(source, lib / name)
            bundled.append(name)

    # Debian/Ubuntu's GCC copyright file contains the full Runtime Library
    # Exception. Copy separate GPL text as well, since the file references it.
    override = Path(os.environ.get("ALTANN_RUNTIME_LICENSE_DIR", "/nonexistent"))
    gcc_copyright = copy_license(
        [override / "libgomp-copyright", override / "copyright",
         Path("/usr/share/doc/libgomp1/copyright")],
        licenses / "libgomp-copyright",
        "Install the libgomp1 documentation, or set ALTANN_RUNTIME_LICENSE_DIR "
        "to a directory containing libgomp-copyright and GPL-3.",
    )
    copyright_text = gcc_copyright.read_text(errors="replace")
    if "GCC RUNTIME LIBRARY EXCEPTION" not in copyright_text:
        copy_license([override / "COPYING.RUNTIME"], licenses / "COPYING.RUNTIME",
                     "Provide COPYING.RUNTIME in ALTANN_RUNTIME_LICENSE_DIR.")
    copy_license([override / "GPL-3", override / "COPYING3",
                  Path("/usr/share/common-licenses/GPL-3")], licenses / "GPL-3",
                 "Install the common-licenses package or provide GPL-3 in ALTANN_RUNTIME_LICENSE_DIR.")
    for soname, package in [("libstdc++.so.6", "libstdc++6"), ("libgcc_s.so.1", "libgcc-s1")]:
        if soname in bundled:
            copy_license([override / f"{package}-copyright",
                          Path("/usr/share/doc") / package / "copyright", gcc_copyright],
                         licenses / f"{package}-copyright",
                         f"Provide the {package} copyright notice in ALTANN_RUNTIME_LICENSE_DIR.")

    resolved = linux_dependencies(core)
    for name in bundled:
        if resolved.get(name, Path("/missing")).resolve() != (lib / name).resolve():
            raise RuntimeError(f"{name} does not resolve inside the bundle. Install altann-core "
                               "with CMake INSTALL_RPATH='$ORIGIN/../lib', then package again.")
    return bundled


def mac_dependencies(binary):
    """Return LC_LOAD_DYLIB names (otool includes a dylib's own ID separately)."""
    dependencies = []
    for line in run("otool", "-L", binary).splitlines()[1:]:
        name = line.strip().split(" (compatibility version", 1)[0]
        if name:
            dependencies.append(name)
    return dependencies


def mac_rpaths(binary):
    text = run("otool", "-l", binary)
    return re.findall(r"cmd LC_RPATH\s+cmdsize \d+\s+path (.*?) \(offset", text)


def brew_libomp_prefix():
    if shutil.which("brew"):
        try:
            prefix = Path(run("brew", "--prefix", "libomp").strip())
            if prefix.is_dir():
                return prefix
        except RuntimeError:
            pass
    for prefix in (Path("/opt/homebrew/opt/libomp"), Path("/usr/local/opt/libomp")):
        if prefix.is_dir():
            return prefix
    return None


def find_libomp(core, dependency, brew_prefix):
    candidates = []
    if dependency.startswith("/"):
        candidates.append(Path(dependency))
    elif dependency.startswith("@loader_path/"):
        candidates.append(core.parent / dependency.removeprefix("@loader_path/"))
    elif dependency.startswith("@rpath/"):
        for rpath in mac_rpaths(core):
            expanded = rpath.replace("@loader_path", str(core.parent)).replace("@executable_path", str(core.parent))
            candidates.append(Path(expanded) / dependency.removeprefix("@rpath/"))
    if brew_prefix:
        candidates.append(brew_prefix / "lib" / "libomp.dylib")
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise RuntimeError("Cannot locate libomp.dylib. Install libomp on the build machine with brew install libomp.")


def bundle_macos(tree, licenses):
    """Rewrite OpenMP loading to the bundled dylib and refresh ad-hoc signatures."""
    core = tree / "bin" / "altann-core"
    dependencies = mac_dependencies(core)
    omp_dependencies = [name for name in dependencies if Path(name).name.startswith("libomp")]
    if len(omp_dependencies) != 1:
        raise RuntimeError("Expected exactly one libomp dependency in altann-core; build with Apple Clang and libomp.")
    prefix = brew_libomp_prefix()
    original = find_libomp(core, omp_dependencies[0], prefix)
    lib = tree / "lib"
    lib.mkdir(exist_ok=True)
    target = lib / "libomp.dylib"
    if original != target.resolve():
        shutil.copy2(original, target)
    # Homebrew files can be read-only; only the staging copy is made writable.
    target.chmod(target.stat().st_mode | 0o200)
    core.chmod(core.stat().st_mode | 0o200)
    run("install_name_tool", "-change", omp_dependencies[0], "@rpath/libomp.dylib", core)
    run("install_name_tool", "-id", "@rpath/libomp.dylib", target)
    rpaths = mac_rpaths(core)
    if "@loader_path/../lib" not in rpaths:
        raise RuntimeError("altann-core lacks @loader_path/../lib in its installed RPATH. Reinstall with the release CMake configuration.")
    # Homebrew may append an absolute library directory during installation.
    # Remove those fallback search paths from the staged executable, ensuring
    # @rpath always selects the bundled OpenMP library after relocation.
    for rpath in rpaths:
        if rpath != "@loader_path/../lib":
            run("install_name_tool", "-delete_rpath", rpath, core)
    if mac_rpaths(core) != ["@loader_path/../lib"]:
        raise RuntimeError("The staged executable still contains an external runtime search path.")

    roots = [original.parent.parent]
    if prefix:
        roots.extend([prefix, prefix.resolve()])
    if os.environ.get("ALTANN_RUNTIME_LICENSE_DIR"):
        roots.insert(0, Path(os.environ["ALTANN_RUNTIME_LICENSE_DIR"]))
    candidates = [root / relative for root in roots for relative in
                  ("LICENSE.txt", "LICENSE.TXT", "LICENSE", "libomp-LICENSE.txt",
                   "share/doc/libomp/LICENSE.txt", "share/libomp/LICENSE.txt")]
    copy_license(candidates, licenses / "libomp-LICENSE.txt",
                 "Provide the installed LLVM OpenMP LICENSE.txt in ALTANN_RUNTIME_LICENSE_DIR; "
                 "the packaging script never downloads licenses.")

    # Reject accidental references to Homebrew/LLVM build-machine libraries.
    # System libc++ and libSystem remain supplied by macOS.
    for binary in (core, target):
        for dependency in mac_dependencies(binary):
            if dependency == "@rpath/libomp.dylib":
                continue
            if not dependency.startswith(("/usr/lib/", "/System/Library/")):
                raise RuntimeError(f"Unbundled macOS dependency in {binary.name}: {dependency}")
    for binary in (target, core):
        run("codesign", "--force", "--sign", "-", binary)
        run("codesign", "--verify", "--strict", binary)
    return ["libomp.dylib"]


def package(prefix, output):
    system = platform.system().lower()
    architecture = {"x86_64": "amd64", "amd64": "amd64", "aarch64": "arm64",
                    "arm64": "arm64"}.get(platform.machine().lower())
    if system not in {"linux", "darwin"} or not architecture:
        raise RuntimeError("Supported package hosts are Linux/macOS on amd64/arm64.")
    label = "macos" if system == "darwin" else "linux"
    archive_root = f"altann-{label}-{architecture}"
    prefix = prefix.resolve()
    for relative in ("bin/altann-core", "bin/altann", "share/altann/altann/__init__.py"):
        if not (prefix / relative).is_file():
            raise RuntimeError(f"Installed file is missing: {prefix / relative}. Run cmake --install first.")
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="altann-package-") as temporary:
        tree = Path(temporary) / archive_root
        shutil.copytree(prefix, tree)
        for name in ("README.md", "LICENSE", "THIRD_PARTY.md", "AUTHORS.md", "CHANGELOG.md"):
            shutil.copy2(ROOT / name, tree / name)
        shutil.copytree(ROOT / "docs", tree / "docs", dirs_exist_ok=True)
        licenses = tree / "licenses"
        licenses.mkdir(exist_ok=True)
        runtimes = (bundle_linux if system == "linux" else bundle_macos)(tree, licenses)
        run(tree / "bin" / "altann-core", "--help", env=clean_runtime_environment())
        (tree / "BUNDLED_RUNTIMES.json").write_text(json.dumps(
            {"platform": label, "architecture": architecture, "libraries": runtimes,
             "licenses": "licenses/", "runtime_requirements": ["Python 3.9+"]},
            indent=2) + "\n")
        # Publish only a complete archive; preserve an existing release on error.
        with tempfile.NamedTemporaryFile(prefix=output.name + ".", suffix=".tmp",
                                         dir=output.parent, delete=False) as handle:
            staged_archive = Path(handle.name)
        try:
            with tarfile.open(staged_archive, "w:gz", dereference=True) as archive:
                archive.add(tree, arcname=archive_root)
            staged_archive.replace(output)
        finally:
            staged_archive.unlink(missing_ok=True)
    print(f"Created {output} (root: {archive_root}; bundled: {', '.join(runtimes)})")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", required=True, type=Path, help="Installed AltAnn tree from cmake --install")
    parser.add_argument("--output", required=True, type=Path, help="Destination .tar.gz archive")
    arguments = parser.parse_args()
    try:
        package(arguments.prefix, arguments.output)
    except (RuntimeError, OSError) as error:
        parser.exit(1, f"package_release: {error}\n")


if __name__ == "__main__":
    main()
