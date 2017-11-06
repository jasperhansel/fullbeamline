


all: pgbind.exe libcpgplot.lib  libcpgplot_g.lib clean

pgbind.exe : pgbind.c
	@cl pgbind.c

libcpgplot.lib:
	@./pgbind ms -h -w pgbind_prototypes
	@for f in `ls cpg*.c`;  do cl /c $$f; done
	@lib /out:libcpgplot.lib cpg*.obj
	@rm cpg*.obj

libcpgplot_g.lib: libcpgplot.lib
	@for f in `ls cpg*.c`;  do cl /c /Zi $$f; done
	@lib /out:libcpgplot_g.lib cpg*.obj

clean:
	@rm cpg*.obj
	@rm cpg*.c
	@rm pgbind.exe

PHONY: clean