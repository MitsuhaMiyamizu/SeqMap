seqmap: 
	g++ -O3 -m64 -o seqmap match.cpp
collect_output:
	g++ -O3 -m64 -o collect_output collect_output.cpp

clean:
	rm -f seqmap 
	rm -f collect_output
