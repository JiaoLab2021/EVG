name: CI

on:
  push:
    branches:
      - main
  pull_request:

env:
  BUILD_TYPE: Release

jobs:
  build:
    runs-on: ubuntu-latest
    
    env:
      DATA_DIR: ./test

    strategy:
      matrix:
        cmake-version: [3.12]
        compiler: [gcc-9]

    steps:
    - uses: actions/checkout@v3
      with:
        submodules: 'recursive'
      env:
        NODE_VERSION: 16.x

    - name: Install dependencies
      run: |
        sudo apt-get update
        sudo apt-get install -y libz-dev

    - name: Build project
      run: |
        cmake . -DCMAKE_BUILD_TYPE=${{ env.BUILD_TYPE }}
        make

    - name: Run tests
      run: |
        echo "Testing"
        ./graphvcf count -v $DATA_DIR/test.vcf.gz
        ./fastAQ count -i $DATA_DIR/test.fa
        python EVG.py -h
