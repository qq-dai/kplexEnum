# The Soruce Code for Enumerating All Maximal k-plexes. 

Compile: 
```
cd /src
make 
```
Note: this code requires g++ version 5.4 or higher to compile.

Example: 
```
./allplexes datas/as-caida2007.txt -k=2 -q=10 
```

- "-k=" is to set the k value for k-plex. 
- "-q=" is a size-constraint of the k-plexes to be enumerated.
