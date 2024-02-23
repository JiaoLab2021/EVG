#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import subprocess
import tempfile
import os
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
