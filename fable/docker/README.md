
docker to run fable FORTRAN converter

# build fable docker
```
docker buildx build -t fable-dock .
```

# run bash in fable docker
# overriding entrypoint
```
docker run -i --entrypoint /usr/bin/env -t fable-dock bash
```

# run fable, see help
```
docker run -it fable-dock
```

# run fable --example to current dir
# fable.cout and fable_cout.cpp are created
```
docker run -i -w $(pwd) --volume=$(pwd):$(pwd):rw -t fable-dock --example
```

# run fable on local files
# send cout to cpp file
```
declare -ar cmda=(
    docker run -i -w $(pwd)
    --volume=$(pwd):$(pwd):rw
    -t fable-dock
     <user fortran file>
     <user fortran file>
    --namespace <user namespace>
)

${cmda[@]} > <user cpp file>
```
