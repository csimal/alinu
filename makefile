# Makefile pour le projet
# Pour compiler, il suffit d'ecrire 'make' dans la shell
# make va compiler tous les fichiers necessaires et les linker dans l'executable
# seuls les fichiers modifies sont recompiles

# le compilateur a utiliser
CC=gcc
# les options du compilateur e.g. -Wall
CFLAGS= -Wall #-g -O0
# les librairies a linker
LIBS= -lm -llapacke

# les objets
OBJS=matrix.o conjugate_gradient.o main.o

# Regle implicite pour compiler les fichiers .c en fichier .o (objets)
.c.o:
	$(CC) $(CFLAGS) -c $< $(LIBS)

main: $(OBJS)
	$(CC) -o main $(OBJS) $(LIBS)

clean:
	rm -f *.o

# Do not delete this line
