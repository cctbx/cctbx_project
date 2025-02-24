
docker to run fable FORTRAN converter

# remember name of docker to build
```
export DOCK_NAME='fable-dock'
```

# build fable docker
```
docker buildx build -t ${DOCK_NAME} .
```

# run bash in fable docker
# overriding entrypoint
```
docker run -i --entrypoint /usr/bin/env -t ${DOCK_NAME} bash
```

# run fable, see help
```
docker run -it ${DOCK_NAME}
```

# run fable --example to current dir
# fable.cout and fable_cout.cpp are created
```
docker run -i -w $(pwd) --volume=$(pwd):$(pwd):rw -t ${DOCK_NAME} --example
```

# run fable on local files
# send cout to cpp file
```
declare -ar cmda=(
    docker run -i -w $(pwd)
    --volume=$(pwd):$(pwd):rw
    -t ${DOCK_NAME}
     <user fortran file>
     <user fortran file>
    --namespace <user namespace>
)

${cmda[@]} > <user cpp file>
```
