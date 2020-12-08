# Here we submit a group of jobs where we fix parameter2 and parameter3
# but want to run with a grid of parameter1 values (1 .. 100)
# also we want to run the same job for all the autosomes

 
param1s=`seq 1 100`
param2=12
param3=0.1

for chr in `seq 1 22`
do
  for param1 in $param1s
  do
    qsub -v PARAMETER1=$param1,PARAMETER2=$param2,PARAMETER3=$param3,CHR_NUM=$chr run.pbs
  done
done

