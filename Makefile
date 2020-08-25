compile:
	g++ -o dc/no_change/dc_no_change dc/no_change/src.cpp

run_dc_no_change:
	g++ -o dc/no_change/dc_no_change dc/no_change/src.cpp
	./dc/no_change/dc_no_change

run_dc_change:
	g++ -o dc/change/dc_change dc/change/src.cpp
	rm -f dc/change/bin/* dc/change/output/test_c/*
	./dc/change/dc_change
	python3 read_test_c.py

clean:
	rm -f dc/no_change/dc_no_change
	rm -f dc/change/dc_change
