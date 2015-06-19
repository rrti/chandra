CC = g++
CFLAGS = -Wall -Wextra

TARGET_REL_ST = chandra-rel-st
TARGET_REL_MT = chandra-rel-mt
TARGET_DBG_ST = chandra-dbg-st
TARGET_DBG_MT = chandra-dbg-mt



chandra-rel-st: main-rel-st.o  mcrtsim-rel-st.o  slab-rel-st.o  grid-rel-st.o  photon-rel-st.o tdf-stackparser-rel-st.o
	$(CC) $(CFLAGS) -O2 -o $(TARGET_REL_ST) *.o -lSDL

chandra-rel-mt: main-rel-mt.o  mcrtsim-rel-mt.o  slab-rel-mt.o  grid-rel-mt.o  photon-rel-mt.o tdf-stackparser-rel-mt.o
	$(CC) $(CFLAGS) -O2 -o $(TARGET_REL_MT) *.o -lSDL -lboost_thread-mt

chandra-dbg-st: main-dbg-st.o  mcrtsim-dbg-st.o  slab-dbg-st.o  grid-dbg-st.o  photon-dbg-st.o tdf-stackparser-dbg-st.o
	$(CC) $(CFLAGS) -o $(TARGET_DBG_ST) *.o -lSDL

chandra-dbg-mt: main-dbg-mt.o  mcrtsim-dbg-mt.o  slab-dbg-mt.o  grid-dbg-mt.o  photon-dbg-mt.o tdf-stackparser-dbg-mt.o
	$(CC) $(CFLAGS) -o $(TARGET_DBG_MT) *.o -lSDL -lboost_thread-mt



## single-threaded release build
main-rel-st.o: Main.cpp
	$(CC) $(CFLAGS) -O2 -c $<
mcrtsim-rel-st.o: MCRTSim.cpp
	$(CC) $(CFLAGS) -O2 -c $<
slab-rel-st.o: Slab.cpp
	$(CC) $(CFLAGS) -O2 -c $<
grid-rel-st.o: Grid.cpp
	$(CC) $(CFLAGS) -O2 -c $<
photon-rel-st.o: Photon.cpp
	$(CC) $(CFLAGS) -O2 -c $<
tdf-stackparser-rel-st.o: TDFStackParser.cpp
	$(CC) $(CFLAGS) -O2 -c $<

## multi-threaded release build
main-rel-mt.o: Main.cpp
	$(CC) $(CFLAGS) -O2 -c $< -DTHREADED
mcrtsim-rel-mt.o: MCRTSim.cpp
	$(CC) $(CFLAGS) -O2 -c $< -DTHREADED
slab-rel-mt.o: Slab.cpp
	$(CC) $(CFLAGS) -O2 -c $< -DTHREADED
grid-rel-mt.o: Grid.cpp
	$(CC) $(CFLAGS) -O2 -c $< -DTHREADED
photon-rel-mt.o: Photon.cpp
	$(CC) $(CFLAGS) -O2 -c $< -DTHREADED
tdf-stackparser-rel-mt.o: TDFStackParser.cpp
	$(CC) $(CFLAGS) -O2 -c $< -DTHREADED

## single-threaded debug build
main-dbg-st.o: Main.cpp
	$(CC) $(CFLAGS) -g -c $<  -DDEBUG
mcrtsim-dbg-st.o: MCRTSim.cpp
	$(CC) $(CFLAGS) -g -c $<  -DDEBUG
slab-dbg-st.o: Slab.cpp
	$(CC) $(CFLAGS) -g -c $<  -DDEBUG
grid-dbg-st.o: Grid.cpp
	$(CC) $(CFLAGS) -g -c $<  -DDEBUG
photon-dbg-st.o: Photon.cpp
	$(CC) $(CFLAGS) -g -c $<  -DDEBUG
tdf-stackparser-dbg-st.o: TDFStackParser.cpp
	$(CC) $(CFLAGS) -g -c $<  -DDEBUG

## multi-threaded debug build
main-dbg-mt.o: Main.cpp
	$(CC) $(CFLAGS) -g -c $<  -DTHREADED -DDEBUG
mcrtsim-dbg-mt.o: MCRTSim.cpp
	$(CC) $(CFLAGS) -g -c $<  -DTHREADED -DDEBUG
slab-dbg-mt.o: Slab.cpp
	$(CC) $(CFLAGS) -g -c $<  -DTHREADED -DDEBUG
grid-dbg-mt.o: Grid.cpp
	$(CC) $(CFLAGS) -g -c $<  -DTHREADED -DDEBUG
photon-dbg-mt.o: Photon.cpp
	$(CC) $(CFLAGS) -g -c $<  -DTHREADED -DDEBUG
tdf-stackparser-dbg-mt.o: TDFStackParser.cpp
	$(CC) $(CFLAGS) -g -c $<  -DTHREADED -DDEBUG



clean:
	rm -f *.o *.dat

