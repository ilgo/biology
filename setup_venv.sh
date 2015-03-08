#!/bin/bash

ENV='venv_bio'

mkdir $ENV
pyvenv-3.4 --without-pip $ENV
. $ENV/bin/activate
curl https://bootstrap.pypa.io/get-pip.py | python

deactivate
. $ENV/bin/activate

pip install numpy==1.8.1
pip install -r requirements.txt
