.PHONY: run reset
run: 
	touch fGenerator5.txt
	gcc fulminePlusPlus.cpp -lstdc++ 
	./a.out
	python KG2.py
reset:  
	rm fGenerator5.txt