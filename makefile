all: 
	make clean laplace
#	make clean yukawa

laplace:
	cd ./src; make lib_lap
	cp ./src/build/libadap_laplace.a ./lib/
	cd ./test; make clean lap laptest3 #laptest6

#yukawa:
#	cd ./src; make lib_yuk
#	cp ./src/build/libadap_yukawa.a ./lib/
#	cd ./test; make clean yuk yuktest3 yuktest6

clean:
	rm -f *~
	(cd ./src; make clean)
	(cd ./example; make clean)
	(cd ./test; make clean)
	(cd ./lib; rm -f *.a)

testall:
	make all | grep "\*\*\* "

profilelap:
	make clean
	cd ./src; make lib_lap PROFILE=-pg
	cp ./src/build/libadap_laplace.a ./lib/
	cd ./test; make clean lap PROFILE=-pg
	rm -f gmon.out
	./test/test_adap_fmm -n 3000000 -d 3 -a 3 -s 30

#profileyuk:
#	make clean
#	cd ./src; make lib_yuk PROFILE=-pg
#	cp ./src/build/libadap_yukawa.a ./lib/
#	cd ./test; make clean yuk PROFILE=-pg
#	rm -f gmon.out
#	./test/test_adap_fmm -b 0.1 -n 3000000 -d 3 -a 3 -s 30


