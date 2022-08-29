#!/bin/bash

necessary_images=("rnaseq")

for image in ${necessary_images[@]};
do
    if [[ "$(docker images -q ${image} 2> /dev/null)" == "" ]];
    then
        echo "### Building image: ${image}:latest ###"
        docker build -t ${image}:latest -f ${image}.Dockerfile .
    fi
done

docker builder prune -f