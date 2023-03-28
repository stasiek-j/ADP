#!/bin/bash

python scripts/download.py data/ -k -p -s -q scripts/organisms.json
python scripts/plots.py data/ --sa --sc --sd --fa --fc --fbc --save_to=img/
python scripts/tables.py data/ --fa --fnt --fc  --sb --format=latex --save_to=tables/

