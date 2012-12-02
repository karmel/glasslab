#!/bin/bash

if [ "$CURRENT_PATH" == "" ]; then
	export CURRENT_PATH="/Users/karmel/Desktop/Projects/GlassLab/Repositories/glasslab/glasslab"
fi 

export PYTHONPATH=$CURRENT_PATH:$PYTHONPATH

export DJANGO_SETTINGS_MODULE=glasslab.config.django_settings

python $CURRENT_PATH/glasslab/inbredstrains/compete/compete_tf_graphs.py $@