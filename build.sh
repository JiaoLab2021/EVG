#!/bin/bash

cmake -DCMAKE_INSTALL_PREFIX=$PREFIX .
make install -j${CPU_COUNT}
