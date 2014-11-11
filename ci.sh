#!/bin/bash

export PYTHONPATH=`pwd`:$PYTHONPATH

# Set up virtual environment
[ ! -e env ] || rm -rf env
virtualenv --python=python2.7 env
. env/bin/activate
# Install testing packages
pip install pep8 pylint pytest pytest-cov

# Install requirements
pip install numpy gr glfw PyOpenGL


# PEP8 conformance check (will cause build failure if python files are not PEP8 conformant)
pep8 $(find . -path ./env -prune -o -name "*.py" -print) | tee pep8.log
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
	echo "Failed PEP8 conformance check!";
	exit 1;
fi


# PyLint code style check (will NOT cause build failure as pylint tends to overreact)
if [[ ! -f pylintrc ]]; then
	echo "pylintrc missing!";
	exit 1;
fi
pylint --rcfile=pylintrc "--msg-template={path}:{line}: [{msg_id}({symbol}), {obj}] {msg}" $(find . -path ./env -prune -o -path ./tests -prune -o -name "*.py" -print) | tee pylint.log


# Run py.test (will cause build failure if a test fails)
if [[ ! -f pytest.ini ]]; then
        echo "pytest.ini missing!";
        exit 1;
fi
if [[ ! -f coveragerc ]]; then
        echo "coveragerc missing!";
        exit 1;
fi
py.test --cov-config coveragerc --cov-report html --cov . | tee py.test.log
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        echo "Failed tests!";
        exit 1;
fi


# Report build success
exit 0
