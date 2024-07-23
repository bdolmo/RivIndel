import os
import subprocess
import sys

# Define the directories for the source packages
HTSLIB_DIR = os.path.join(os.getcwd(), "bcftools-1.20/htslib-1.20")
BCFTOOLS_DIR = os.path.join(os.getcwd(), "bcftools-1.20")

def run_command(command, cwd=None, env=None):
    """Run a shell command."""
    print(f"Running: {' '.join(command)}")
    subprocess.check_call(command, cwd=cwd, env=env)


def download_and_extract(url, dest):
    """Download and extract a tarball."""
    tarball = os.path.join(dest, os.path.basename(url))
    run_command(["wget", url, "-O", tarball])
    if tarball.endswith(".tar.gz"):
        run_command(["tar", "-xzf", tarball, "-C", dest])
    elif tarball.endswith(".tar.bz2"):
        run_command(["tar", "-xjf", tarball, "-C", dest])


def install_htslib(prefix):
    """Install HTSlib from the provided source directory."""
    env = os.environ.copy()
    env['CFLAGS'] = f"-I{prefix}/include"
    env['LDFLAGS'] = f"-L{prefix}/lib"
    run_command(["autoheader"], cwd=HTSLIB_DIR)
    run_command(["autoconf"], cwd=HTSLIB_DIR)
    run_command(["./configure", "--disable-shared", f"--prefix={prefix}"], cwd=HTSLIB_DIR, env=env)
    run_command(["make"], cwd=HTSLIB_DIR)
    run_command(["make", "install"], cwd=HTSLIB_DIR)

def install_bcftools(prefix):
    """Install BCFtools from the provided source directory."""
    env = os.environ.copy()
    env['CFLAGS'] = f"-I{prefix}/include"
    env['LDFLAGS'] = f"-L{prefix}/lib"
    run_command(["autoheader"], cwd=BCFTOOLS_DIR)
    run_command(["autoconf"], cwd=BCFTOOLS_DIR)
    run_command(["./configure", f"--with-htslib={HTSLIB_DIR}", "--disable-plugins", f"--prefix={prefix}"], cwd=BCFTOOLS_DIR, env=env)
    run_command(["make"], cwd=BCFTOOLS_DIR)
    run_command(["make", "install"], cwd=BCFTOOLS_DIR)

def main():
    install_dir = os.path.join(os.getcwd(), "install_dir")
    os.makedirs(install_dir, exist_ok=True)

    try:
        install_htslib(install_dir)
        install_bcftools(install_dir)
        print("Installation of bcftools and dependencies completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
