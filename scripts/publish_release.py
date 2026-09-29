"""Publish one tested main-branch version, without overwriting a public release.

This runs only after all four native CI jobs succeed. Version changes are
reviewed in ordinary PRs; the workflow creates the version tag automatically.
The GitHub CLI reads GH_TOKEN from its environment and never prints the token.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess


def run(*args, check=True):
    return subprocess.run(args, text=True, capture_output=True, check=check)


def main():
    root = Path(__file__).resolve().parents[1]
    version = re.search(r"^__version__ = ['\"](\d+\.\d+\.\d+)['\"]$",
                        (root / "altann/__init__.py").read_text(), re.MULTILINE)
    if not version:
        raise RuntimeError("Expected a three-component version in altann/__init__.py")
    version = version.group(1)
    if f"project(AltAnn VERSION {version} " not in (root / "CMakeLists.txt").read_text():
        raise RuntimeError("CMake and Python versions disagree")
    tag = f"v{version}"
    sha = os.environ["GITHUB_SHA"]
    repo = os.environ["GH_REPO"]
    # A published release is immutable to this automation. Main can receive
    # documentation/maintenance changes without generating another version.
    found = run("gh", "release", "view", tag, "--json", "isDraft", check=False)
    if found.returncode == 0 and not json.loads(found.stdout)["isDraft"]:
        print(f"{tag} is already published; no release changes")
        return
    existing_tag = run("git", "rev-parse", "--verify", f"refs/tags/{tag}^{{commit}}", check=False)
    if existing_tag.returncode == 0 and existing_tag.stdout.strip() != sha:
        raise RuntimeError(f"{tag} identifies a different commit; bump the version instead")

    artifacts = [root / "dist" / f"altann-{os_name}-{arch}.tar.gz"
                 for os_name in ("linux", "macos") for arch in ("amd64", "arm64")]
    for artifact in artifacts:
        if not artifact.is_file() or not artifact.stat().st_size:
            raise RuntimeError(f"Missing release artifact: {artifact}")
    checksum = root / "dist" / "SHA256SUMS"
    checksum.write_text("".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n"
                               for p in artifacts))
    if existing_tag.returncode != 0:
        # This lightweight tag points to the exact commit whose artifacts passed
        # CI. We publish below in the same job; no tag-triggered workflow is needed.
        run("gh", "api", f"repos/{repo}/git/refs", "--method", "POST",
            "-f", f"ref=refs/tags/{tag}", "-f", f"sha={sha}")
    if found.returncode != 0:
        run("gh", "release", "create", tag, "--verify-tag", "--draft",
            "--title", f"AltAnn {tag}", "--generate-notes")
    run("gh", "release", "upload", tag, *(str(p) for p in [*artifacts, checksum]), "--clobber")
    run("gh", "release", "edit", tag, "--draft=false")
    print(f"Published {tag} with all four platform packages and SHA256SUMS")


if __name__ == "__main__":
    main()
