compile:
	g++ -o dc_fix/exec dc_fix/src.cpp

run_dc_fix:
	g++ -o dc_fix/exec dc_fix/src.cpp
	./dc_fix/exec

run_dc_dev:
	g++ -o dc_dev/exec dc_dev/src.cpp
	rm -f dc_dev/bin/* dc_dev/output/conc/*
	./dc_dev/exec
	python3 read_test_c.py

clean:
	rm -f dc_fix/exec
	rm -f dc_dev/exec
