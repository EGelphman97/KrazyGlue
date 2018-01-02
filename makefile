.PHONY: run reset
run: 
	touch fGenerator5.txt
	touch outputH.txt
	gcc fulminePlusPlus.cpp -lstdc++ 
	./a.out
	python KG2.py
reset:  
	rm fGenerator5.txt
	rm outputH.txt