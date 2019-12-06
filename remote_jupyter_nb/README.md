# Workflow

This directory is more about how to open jupyter notebook in a remote server and use the local browser to skim it.

The tips and steps are from this nice [post](https://ljvmiranda921.github.io/notebook/2018/01/31/running-a-jupyter-notebook/)

## On the remote side

```
$ jupyter notebook --no-browser --port=XXXX
```

## On the local machine 

```
$ ssh -N -f -L localhost:YYYY:localhost:XXXX remoteuser@remotehost
```

# An working example:

On remote side

```
$ jupyter notebook --no-browser --port=8888
```

On local machine

```
$ ssh -N -f -L localhost:YYYY:localhost:XXXX yanyul@server
```