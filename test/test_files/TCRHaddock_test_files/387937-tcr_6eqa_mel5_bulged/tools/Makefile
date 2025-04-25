#####################################
#                                   #
#  HADDOCK2.0 tools makefile        #
#                                   #
#####################################

CPP=g++
CPFLAGS=-O2

EXEC=cluster_struc contact contact-chainID haddock-decompress-fastfunc contact_fcc contact_fcc_lig

all: 
	make $(EXEC)

cluster_struc: cluster_struc.cpp
	$(CPP) $(CPFLAGS) -o cluster_struc cluster_struc.cpp
	   
contact: contact.cpp
	$(CPP) $(CPFLAGS) -o contact contact.cpp

contact-chainID: contact-chainID.cpp
	$(CPP) $(CPFLAGS) -o contact-chainID contact-chainID.cpp

contact_fcc: contact_fcc.cpp
	$(CPP) $(CPFLAGS) -o contact_fcc contact_fcc.cpp

contact_fcc_lig: contact_fcc_lig.cpp
	$(CPP) $(CPFLAGS) -o contact_fcc_lig contact_fcc_lig.cpp
	   
haddock-decompress-fastfunc: haddock-decompress-fastfunc.cpp
	$(CPP) $(CPFLAGS) -o haddock-decompress-fastfunc haddock-decompress-fastfunc.cpp
	   
clean :
	/bin/rm -f $(EXEC)

