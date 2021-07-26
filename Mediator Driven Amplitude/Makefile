#make file for mediator dependent tree amplitude

run.exe: med_hel_amp.o Spinor_functions.o functions_file.o
	g++ med_hel_amp.o Spinor_functions.o functions_file.o -o run

med_hel_amp.o: med_hel_amp.cpp Spinor.h constants.h structures_list.h functions_list.h
	g++ -c med_hel_amp.cpp

Spinor_functions.o: Spinor_functions.cpp Spinor.h constants.h structures_list.h
	g++ -c Spinor_functions.cpp

functions_file.o: functions_file.cpp
	g++ -c functions_file.cpp

clean:
	rm -f *.o run.exe