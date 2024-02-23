#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import run_cmd


# check file
def check_software_existence(software_name, env_path):
    cmd = f"which {software_name}"
    stdout, stderr, log_out = run_cmd.run(cmd, env_path)
    return stdout, stderr, log_out
