.PHONY: tables quick clean

tables:
	Rscript run_all.R

quick:
	Rscript run_all.R --quick

clean:
	rm -f output/*.tex
