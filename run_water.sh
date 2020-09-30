#!/bin/bash

docker run -it -v $PWD/.config:/root/.config scl3/task_hii_water python hii_water.py -r 'Afrotropic'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_water python hii_water.py -r 'Australasia'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_water python hii_water.py -r 'Indomalayan'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_water python hii_water.py -r 'Nearctic'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_water python hii_water.py -r 'Neotropic'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_water python hii_water.py -r 'Oceania'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_water python hii_water.py -r 'Palearctic'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_water python hii_water.py -r 'HighArctic'
