"""Shared pieces of the EVG pipeline: program paths, command runner with resume records, and job pool."""

import concurrent.futures
import glob
import hashlib
import json
import logging
import os
import re
import subprocess
import threading
import uuid

# Directory holding EVG.py, graphvcf, fastAQ, src/evg_pipeline/ and include/ (same in the source tree and in bin/).
CODE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# External commands see only PATH.
ENV = {"PATH": os.environ.get("PATH")}
# EVG's own binaries are not recorded: rebuilding EVG does not re-run finished commands.
OWN_BINARIES = (os.path.join(CODE_DIR, "graphvcf"), os.path.join(CODE_DIR, "fastAQ"))
# Inside a run_jobs task, JOB.last is the record id of the job's previous command.
JOB = threading.local()
# Set after the first failed job: running jobs stop before their next command.
STOP = threading.Event()

logger = logging.getLogger("EVG")


def stat_files(cmd, cwd, files):
    """Return {path: [size, mtime_ns]} for the regular files named in cmd (also after "=") or listed in `files` (literal or glob)."""
    paths = []
    for word in re.split(r"[\s'\"`<>|&;()]+", cmd):
        paths += [word, word.split("=", 1)[1]] if "=" in word else [word]
    paths = [os.path.join(cwd, path) for path in paths if path]
    for entry in files:
        paths += sorted(glob.glob(os.path.join(cwd, entry))) if glob.has_magic(entry) else [os.path.join(cwd, entry)]
    stats = {}
    for path in paths:
        path = os.path.normpath(path)
        if path.startswith(("/dev/", "/proc/")) or path in OWN_BINARIES or not os.path.isfile(path):
            continue
        status = os.stat(path)
        stats[path] = [status.st_size, status.st_mtime_ns]
    return stats


def find_change(record, last):
    """Return why a finished command must run again (a changed or missing file, or a re-run previous command), or ""."""
    for path, stat in record["files"].items():
        try:
            status = os.stat(path)
        except OSError:
            return path
        if [status.st_size, status.st_mtime_ns] != stat:
            return path
    if record["after"] != last:
        return "previous command of this job"
    return ""


def run_cmd(cmd, cwd, files=()):
    """Run a shell command in cwd unless it already finished with unchanged files; return its stdout; raise RuntimeError on failure or a missing literal `files` entry."""
    key = re.sub(r"(?<!\S)(-t|-@|-p|-j|--threads)([ =])\d+(?=\s|$)", r"\1\2N", cmd)
    record_file = os.path.join(cwd, ".evg", hashlib.sha256(key.encode()).hexdigest() + ".json")
    last = getattr(JOB, "last", None)

    # ------ Finished before ------
    try:
        with open(record_file) as f:
            record = json.load(f)
    except (OSError, ValueError):
        record = None
    if not isinstance(record, dict) or not {"cmd", "stdout", "files", "id", "after"} <= record.keys():
        record = None
    if record is not None:
        change = find_change(record, last)
        if not change:
            logger.info(f"  - Skipped: {cmd}", stacklevel=2)
            if hasattr(JOB, "last"):
                JOB.last = record["id"]
            return record["stdout"]
        logger.info(f"  - Changed: {change}", stacklevel=2)
        try:
            os.remove(record_file)
        except FileNotFoundError:
            pass

    # ------ Run ------
    if STOP.is_set():
        raise RuntimeError(f"Stopped after an earlier failure: {cmd}")
    logger.info(f"  $ {cmd}", stacklevel=2)
    result = subprocess.run("set -o pipefail; " + cmd, shell=True, executable="/bin/bash", cwd=cwd, env=ENV,
                            capture_output=True)
    stdout = result.stdout.strip().decode("utf-8", errors="replace")
    stderr = result.stderr.strip().decode("utf-8", errors="replace")
    if result.returncode != 0 or "FileNotFoundError" in stderr or "command not found" in stderr \
            or "error" in stderr or "Error" in stderr:
        raise RuntimeError(f"Command failed (exit {result.returncode}): {cmd}\n{stderr}")

    # ------ Record ------
    missing = [entry for entry in files if not glob.has_magic(entry) and not os.path.isfile(os.path.join(cwd, entry))]
    if missing:
        raise RuntimeError(f"Command finished without writing {os.path.join(cwd, missing[0])}: {cmd}")
    record = {"cmd": cmd, "stdout": stdout, "files": stat_files(cmd, cwd, files), "id": uuid.uuid4().hex, "after": last}
    os.makedirs(os.path.dirname(record_file), exist_ok=True)
    temp_file = f"{record_file}.{uuid.uuid4().hex}.tmp"
    with open(temp_file, "w") as f:
        json.dump(record, f)
    os.replace(temp_file, record_file)
    if hasattr(JOB, "last"):
        JOB.last = record["id"]
    return stdout


def write_list(path, text):
    """Write a sample list unless it already holds this text (an unchanged list keeps its mtime for --resume)."""
    if os.path.isfile(path):
        with open(path) as f:
            if f.read() == text:
                return
    with open(path, "w") as f:
        f.write(text)


def get_version(software, cwd):
    """Return the version printed by running the software without arguments in cwd, or "" if none can be read."""
    try:
        result = subprocess.run(software, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                env=ENV, universal_newlines=True, timeout=10)
    except Exception:
        return ""
    match = re.search(r"version[:]?[\s]+(v?[\d+\.]+)", result.stdout + result.stderr, re.IGNORECASE)
    return match.group(1).lstrip("vV") if match else ""


def parse_version(version):
    """Turn a version such as "1.63" into (1, 63, 0) for numeric comparison; text that is not a version gives ()."""
    if not re.match(r"\d", version):
        return ()
    numbers = [int(number) for number in re.findall(r"\d+", version)]
    return tuple(numbers + [0] * (3 - len(numbers)))


def run_jobs(tasks, workers, progress=True):
    """Run (function, *args) tasks in a thread pool and return results in task order; on a failure skip queued tasks, wait for running ones and re-raise."""
    def run_task(function, *args):
        if STOP.is_set():
            raise RuntimeError("Stopped after an earlier failure.")
        JOB.last = None
        try:
            return function(*args)
        except BaseException:
            STOP.set()
            raise
        finally:
            del JOB.last

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
    futures = []
    try:
        futures = [executor.submit(run_task, *task) for task in tasks]
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            future.result()
            if progress and len(futures) > 1:
                logger.info(f"  # {done * 100 // len(futures):3d}% ({done:{len(str(len(futures)))}d}/{len(futures)}) done", stacklevel=2)
    except BaseException as error:
        STOP.set()
        executor.shutdown(wait=False, cancel_futures=True)
        if isinstance(error, KeyboardInterrupt):
            raise
        running = sum(future.running() for future in futures)
        if running:
            logger.warning(f"A job failed; waiting for {running} running job(s) to finish ...", stacklevel=2)
        executor.shutdown()
        raise
    executor.shutdown()
    return [future.result() for future in futures]
