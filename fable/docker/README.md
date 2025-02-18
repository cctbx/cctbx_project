
docker to run fable FORTRAN converter

# build
docker buildx build -t fable-dock --build-arg=BUILD_PARENT_DOCKPATH=python:3.9-slim .

# run bash, overriding entrypoint
docker run -i --entrypoint /usr/bin/env -t fable-dock bash

# run fable
distro.sh run
docker run -it fable-dock

# run fable --example to current dir
# fable.cout and fable_cout.cpp are created
docker run -i -w $(pwd) --volume=$(pwd):$(pwd):rw -t fable-dock --example

# run fable on local files
declare -ar cmda=(
    docker run -i -w $(pwd)
    --volume=$(pwd):$(pwd):rw
    -t fable-dock
    fort1.for
    fort2.for
    --namespace calcs
)

${cmda[@]} > calcs.cpp
