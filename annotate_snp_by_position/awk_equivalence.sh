# ARGS1: file to annotate
# ARGS2: lookup table 

myfile=$1
lookup=$2

awk '{
  if(FNR==NR) {
    memory["chr"$2"-"$3]=$0
  } else {
    key=$1"-"$2; 
    if(key in memory){ 
      print memory[key]"\t"$4
    }
  }
}' <(zcat < $myfile) <(zcat < $lookup)