#!/bin/bash

jupyter nbconvert --to notebook --execute --inplace $1 --ExecutePreprocessor.timeout=-1

