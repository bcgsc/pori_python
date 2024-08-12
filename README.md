
# PORI Python Adaptor

![build](https://github.com/bcgsc/pori_python/workflows/build/badge.svg) [![codecov](https://codecov.io/gh/bcgsc/pori_python/branch/master/graph/badge.svg)](https://codecov.io/gh/bcgsc/pori_python)

This repository is part of the [Platform for Oncogenomic Reporting and Interpretation (PORI)](https://bcgsc.github.io/pori/).

This is a python adaptor package for querying the GraphKB API and IPR API.

This python tool takes in variant inputs as tab-delimited files and annotates them using GraphKB.
The resulting output is uploaded to IPR as a report. Additional report content such as images and
metadata can be passed to be included in the report upload.

For documentation on how to create reports using the IPR adaptor, see the [main documentation site](https://bcgsc.github.io/pori/) for the platform. For the GraphKB adaptor, see the [user manual](https://bcgsc.github.io/pori/graphkb/scripting/).

- [Getting Started](#getting-started)
  - [Install (For developers)](#install-for-developers)
- [Documentation](#documentation)
- [Deployment (Publishing)](#deployment-publishing)

## Getting Started

### Install (For developers)

clone this repository

```bash
git clone https://github.com/bcgsc/pori_python.git
cd pori_python
```

create a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate
```

install the package and its development dependencies

```bash
pip install -U pip setuptools
pip install -e .[dev]
```

Run the tests:

Export usernames, passwords, and set test options.

Note that IPR tests will try to use the BCGSC production GraphKB API by default.
If you want to test interaction with a different instance, you will need to
set the GraphKB variables.

Set EXCLUDE vars to 1 if you don't want to run these tests. 
ONCOKB and BCGSC tests are enabled by default.

```bash
export IPR_USER='pori_admin'
export IPR_PASS='pori_admin'
export IPR_URL='http://localhost:8081/api'
export GRAPHKB_USER='pori_admin'
export GRAPHKB_PASS='pori_admin'
export GRAPHKB_URL='http://localhost:8080/api'
EXCLUDE_BCGSC_TESTS = 1
EXCLUDE_ONCOKB_TESTS = 1
```

If you want to run tests that upload reports to a live IPR instance,
specify the url of the IPR API you want to use and set the test var to 1.
These tests are disabled by default.

```bash
export IPR_TEST_URL='http://localhost:8081/api'
INCLUDE_UPLOAD_TESTS = 1
```

```bash
pytest tests
```

## Documentation

The user documentation for this tool is hosted with the [main documentation site](https://bcgsc.github.io/pori/).

Developers: Any updates to this tool should be edited and reflected in the main site documentation as well.


## Deployment (Publishing)

Install the deployment dependencies

```bash
pip install .[deploy]
```

Build the distribution files

```bash
python setup.py install sdist bdist_wheel
```

Upload the distibutions to the package server (`-r` is defined in your pypirc)

```bash
twine upload -r bcgsc dist/*
```
