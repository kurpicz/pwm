SRCDIR=~/docs/work/gitlab/cocada/src

# pull:
# 	mkdir -p src
# 	rm -f src/*
# 	cp ${SRCDIR}/alphabet.* src/
# 	cp ${SRCDIR}/arrayutil.* src/
# 	cp ${SRCDIR}/binheap.* src/
# 	cp ${SRCDIR}/bitarray.* src/
# 	cp ${SRCDIR}/bitsandbytes.* src/
# 	cp ${SRCDIR}/bitvector.* src/
# 	cp ${SRCDIR}/bytearray.* src/
# 	cp ${SRCDIR}/cocadautil.* src/
# 	cp ${SRCDIR}/csrsbitarray.* src/
# 	cp ${SRCDIR}/cstringutil.* src/
# 	#cp ${SRCDIR}/debruijngraph.* src/
# 	#cp ${SRCDIR}/djb2.* src/
# 	cp ${SRCDIR}/dynarray.* src/
# 	cp ${SRCDIR}/dynstr.* src/
# 	cp ${SRCDIR}/hashmap.* src/
# 	cp ${SRCDIR}/hashset.* src/
# 	cp ${SRCDIR}/huffcode.* src/
# 	#cp ${SRCDIR}/marshall.* src/
# 	cp ${SRCDIR}/mathutil.* src/
# 	#cp ${SRCDIR}/queue.* src/
# 	#cp ${SRCDIR}/srsbitarray.* src/
# 	#cp ${SRCDIR}/stack.* src/
# 	cp ${SRCDIR}/strstream.* src/
# 	#cp ${SRCDIR}/vebwtree.* src/
# 	cp ${SRCDIR}/wavtree.* src/



build: fs1.c fs2.c
	cd src/; rm *.o; gcc -O3 -c *.c; cd ..
	gcc -O3 -I src/ -o fs1 fs1.c src/*.o -lm
	gcc -O3 -I src/ -o fs2 fs2.c src/*.o -lm

install: build

clean:
	rm fs 
