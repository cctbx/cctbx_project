
name=dxtbx

default: $(name).pdf


$(name).pdf: $(name).dvi
	dvipdf $(name).dvi
    
$(name).dvi: $(name).tex
	latex $(name).tex
	bibtex $(name).aux
	latex $(name).tex
	latex $(name).tex
	latex $(name).tex
	
clean:
	rm -vf $(name).aux 
	rm -vf $(name).log 
	rm -vf $(name).dvi 
	rm -vf $(name).pdf 
	rm -vf $(name).out
	rm -vf $(name).bbl
	rm -vf $(name).blg
