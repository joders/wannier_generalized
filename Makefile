
main:	inWan_phony

exec:	inWan_phony

display:	pdf
	./graph.sh

pdf:	output
	./createPlots.sh

#output:	inWan_phony
#	echo > output/dirFiller; rm output/*; ./inWan 1> output/stdout; xmessage job finished; cp t.cpp output/code.cpp
output:	inWan_phony
	if test -d output; then rm -r output; fi; ./inWan 1> /dev/null; xmessage job finished; cp t.cpp output/code.cpp

# always recompile (condition-less)
inWan_phony:	t.cpp streamer.cpp
	g++ -o inWan -larmadillo -lgsl -I include t.cpp

debug:	symbols
	cgdb inWan

symbols:	t.cpp
	/usr/bin/g++ -o inWan -DEXTRACHECKS -Wall -ggdb -larmadillo -lgsl -I include t.cpp 
#				-Wall -Weffc++ -pedantic  \
				-pedantic-errors -Wextra  -Wall -Waggregate-return -Wcast-align \
				-Wcast-qual  -Wchar-subscripts  -Wcomment -Wconversion \
				-Wdisabled-optimization \
				-Werror \
				-Wformat  \
				-Wformat=2 \
				-Wformat-nonliteral -Wformat-security  \
				-Wformat-y2k \
				-Wimport  -Winit-self  -Winline \
				-Winvalid-pch   \
				-Wunsafe-loop-optimizations  -Wlong-long -Wmissing-braces \
				-Wmissing-field-initializers -Wmissing-format-attribute   \
				-Wmissing-include-dirs -Wmissing-noreturn \
				-Wpacked  -Wpadded -Wparentheses  -Wpointer-arith \
				-Wredundant-decls -Wreturn-type \
				-Wsequence-point  -Wshadow -Wsign-compare  -Wstack-protector \
				-Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default \
				-Wswitch-enum -Wtrigraphs  -Wuninitialized \
				-Wunknown-pragmas  -Wunreachable-code -Wunused \
				-Wunused-function  -Wunused-label  -Wunused-parameter \
				-Wunused-value  -Wunused-variable  -Wvariadic-macros \
				-Wvolatile-register-var  -Wwrite-strings \
				-Wfloat-equal

debug_nv:	symbols_nv
	cgdb inWan

symbols_nv:	t.cpp
	/usr/bin/g++ -o inWan -Wall -ggdb -larmadillo -lgsl -I ../../../../../include t.cpp


.PHONY: debug symbols display pdf output
