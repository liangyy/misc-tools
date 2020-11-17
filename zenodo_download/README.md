Sometimes you may want to download Zenodo datasets on remote server/cluster.
So, it is handy to have a directly downloadable URL.

The way to get it is the following:

```
# Use zenodo REST API
curl https://zenodo.org/api/records/[zenodo_id]
# Grasp the URL from "files" -> "links"
wget [dataset-URL]
```