CFLAG = -framework Accelerate
GSL = -lgsl

gr2sf: gr2sf.cc
	g++ -o gr2sf gr2sf.cc $(GSL) $(CFLAG) 
apply_gr2sf:
	g++ -o apply_gr2sf apply_gr2sf.cc $(CFLAG)
hard_disc_rdf:
	g++ -O3 -o hard_disc_rdf hard_disc_rdf.cc $(GSL) $(CFLAG)
