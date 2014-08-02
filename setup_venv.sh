#!/bin/bash

ENV='venv_bio'

mkdir $ENV
virtualenv -p /usr/bin/python3 $ENV
pip install -r requirements.txt
