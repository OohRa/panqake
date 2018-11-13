CC = gfortran
FLIST = panqake0.f90
EXEC = panqake0

panqake0:
	$(CC) $(FLIST) -o $(EXEC)

clean:
	rm -f $(EXEC)
