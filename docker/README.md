# dyngen manuscript in a Docker container

[dynverse/dyngen_manuscript](https://hub.docker.com/r/dynverse/dyngen_manuscript) contains all necessary packages to reproduce the dyngen manuscript from start to finish.

## Running the container
To run the container, you can use the following command.

```sh
docker run --rm -p 127.0.0.1:8787:8787 -e DISABLE_AUTH=true -v `pwd`:/home/rstudio/dyngen_manuscript -v "$HOME/.cache:/home/rstudio/.cache" dynverse/dyngen_manuscript:1.0
```

Keep this window open, and open up a browser and go to [127.0.0.1:8787](127.0.0.1:8787). 

The command can be dissected as follows.

```sh
docker run \

  # remove container after use
  --rm \
  
  # specify which port rstudio server uses
  -p 127.0.0.1:8787:8787 \
  
  # disable authentication because I'm lazy
  -e DISABLE_AUTH=true \
  
  # mount the current directory at '~/dyngen_manuscript'
  -v `pwd`:/home/rstudio/dyngen_manuscript \
  
  # also mount the cache to save some time
  -v "$HOME/.cache:/home/rstudio/.cache" \
  
  # specify which container to run
  dynverse/dyngen_manuscript:1.0
```

## Update the container

If a newer version of the container has been released, you can update it by running the following command.
```sh
docker pull dynverse/dyngen_manuscript
```


## Building the container

To rebuild this docker container from scratch, run the following command in the `dyngen_manuscript` root directory.

```sh
docker build -t dynverse/dyngen_manuscript:1.0 -f docker/Dockerfile --build-arg GITHUB_PAT=$GITHUB_PAT .
```

GITHUB_PAT should be an environment variable corresponding to the Personal Access Token created by following [this tutorial](https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token).


