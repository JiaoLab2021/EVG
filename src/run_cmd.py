#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import subprocess
import tempfile
import os
import re
import logging


# log
logger = logging.getLogger('run_cmd')
formatter = logging.Formatter('[%(asctime)s] %(message)s')
handler = logging.StreamHandler()  # output to the console
handler.setFormatter(formatter)
logger.addHandler(handler)


# run
def run(command, envPath):
    """
    Executes a command in the shell and captures the standard output, standard error, and log.

    Args:
        command (str): The command to be executed.
        envPath (dict): The environment variables to be set for the command.

    Returns:
        tuple: A tuple containing the standard output (str), standard error (str), and log (str).
    """
    # Create output file
    stdout_file = tempfile.NamedTemporaryFile(delete=False)
    stderr_file = tempfile.NamedTemporaryFile(delete=False)

    # submit task
    proc = subprocess.Popen(command, shell=True, stdout=stdout_file, stderr=stderr_file, env=envPath)

    log = f'CMD: {command}'
    logger.error(log)

    # Reset log is used to judge whether to exit normally
    log = ""

    # Wait for the command to complete
    proc.wait()

    # Read the output file and close the file handle
    with open(stdout_file.name, 'rb') as f:
        stdout_data = f.read()
    os.unlink(stdout_file.name)

    with open(stderr_file.name, 'rb') as f:
        stderr_data = f.read()
    os.unlink(stderr_file.name)

    exit_state = proc.returncode
    stdout = stdout_data.strip().decode('utf-8')
    stderr = stderr_data.strip().decode('utf-8')

    # standard output and error output
    if exit_state != 0 or "FileNotFoundError" in stderr or \
            "command not found" in stderr or \
            "error" in stderr or \
            "Error" in stderr:
        log = f'Error occurred with status code: {exit_state}. Error message: {stderr}.'

    return stdout, stderr, log


def get_version(software_name, envPath, args=""):
    """
    Runs a software command with optional arguments in a given environment and extracts version info.

    Args:
        software_name (str): The software command to run (e.g., 'vg', 'PanGenie').
        envPath (dict): The environment in which to run the command.
        args (str): Command-line arguments as a string, e.g. "-h", "--version", or other flags.

    Returns:
        str: Extracted version or full output if version not found.
    
    This function runs the given software command with the provided environment and arguments. It
    attempts to extract version information from the command's output using a regular expression.
    The function handles common output formats and returns the version number or an error message.
    """
    try:
        full_cmd = f"{software_name} {args}".strip()

        result = subprocess.run(
            full_cmd, 
            shell=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            env=envPath, 
            universal_newlines=True, 
            timeout=10
        )

        text = result.stdout + result.stderr

        # Try to extract version using regex
        match = re.search(r'version[:]?[\s]+(v?[\d+\.]+)', text, re.IGNORECASE)
        if match:
            version_raw = match.group(1)
            version_clean = version_raw.lstrip("vV")
            logger.info(f"{software_name} version: {version_clean}")
            return version_clean
        else:
            logger.warning(f"[Warning] Could not parse version info from {software_name}")
            return "0.0.0"
    except subprocess.CalledProcessError as e:
        return f"Error running {software_name}: {e.output.decode('utf-8')}"
    except Exception as e:
        return f"Exception: {str(e)}"