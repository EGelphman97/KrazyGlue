.PHONY: clean run
clean:
	rm test5.txt
	rm output5.txt
run: 
	touch test5.txt
	touch output5.txt
	python KrazyGlue.py
