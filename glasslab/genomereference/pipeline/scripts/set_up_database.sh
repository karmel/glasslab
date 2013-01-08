#!/bin/bash

export PYTHONPATH=$CURRENT_PATH:$PYTHONPATH

export DJANGO_SETTINGS_MODULE=glasslab.config.current_settings

python $CURRENT_PATH/glasslab/genomereference/pipeline/set_up_database.py $@
