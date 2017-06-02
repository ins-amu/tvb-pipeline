#!/bin/bash

if docker info &> /dev/null
then
    docker="docker "
else
    docker="sudo docker "
fi

echo $docker

if $docker info &> /dev/null
then
    $docker run --rm -it maedoc/tvb-make $@
else
    echo "unable to access Docker, is it running?"
    exit 1
fi
