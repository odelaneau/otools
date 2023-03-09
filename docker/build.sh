#Clean up
rm ../static_bins/OTOOLS_*
rm resources/OTOOLS_*

#Compile mendel
cd ../mendel/
make clean
make -j static_exe
cp bin/OTOOLS_mendel_static ../static_bins/.

#Compile swapalleles
cd ../swapalleles/
make clean
make -j static_exe
cp bin/OTOOLS_swapalleles_static ../static_bins/.

#Compile liftover
cd ../liftover/
make clean
make -j static_exe
cp bin/OTOOLS_liftover_static ../static_bins/.

#Buld docker image
LAB=otools_$(git log -1 --format=%cd --date=short)\_$(git rev-parse --short HEAD)

cd ../docker/
mkdir -p resources
cp ../static_bins/OTOOLS* resources/.

docker build -t $LAB -f Dockerfile .
docker save $LAB | gzip -c > $LAB\.tar.gz
